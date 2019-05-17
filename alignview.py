#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import os
import sys
import errno
import logging
import argparse
import collections
import getreads
assert sys.version_info.major >= 3, 'Python 3 required'

DESCRIPTION = """Take an alignment and print a version with a consensus sequence and the conserved
bases masked."""

# The ascii values that represent a 0 PHRED score.
QUAL_OFFSETS = {'sanger':33, 'solexa':64}
DEFAULT_SEQ_COLUMN = 2
DEFAULT_QUAL_COLUMN = 3

def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  input = parser.add_argument_group('Input')
  input.add_argument('input', metavar='align.fa', type=argparse.FileType('r'), default=sys.stdin,
                      nargs='?',
    help='Aligned input sequences. Can be FASTA or just the plain sequences, one per line.')
  input.add_argument('-f', '--format', choices=list(getreads.FORMATS)+['msa'],
    help='Format of the input. Will be inferred from the filename if not given. "lines" is the most '
         'basic format: Just the sequences, one per line.')
  input.add_argument('-S', '--seq-column', type=int,
    help='For tsv (tab-delimited) input, this column contains the sequence data (1-indexed). '
         'If you specify this manually without specifying --qual-column, it will assume there are '
         'no quality scores. '
         'Default: '+str(DEFAULT_SEQ_COLUMN))
  input.add_argument('-Q', '--qual-column', type=int,
    help='For tsv (tab-delimited) input, this column contains the quality scores (1-indexed). '
         'Default: '+str(DEFAULT_QUAL_COLUMN))
  input.add_argument('-F', '--qual-format', choices=QUAL_OFFSETS.keys(), default='sanger',
    help='FASTQ quality score format. Sanger scores are assumed to begin at \'{}\' ({}). '
         'Default: %(default)s.'.format(QUAL_OFFSETS['sanger'], chr(QUAL_OFFSETS['sanger'])))
  cons = parser.add_argument_group('Consensus calling')
  cons.add_argument('-q', '--qual-thres', type=int, default=25,
    help='Quality score threshold. If quality scores are present, don\'t show bases with lower '
         'quality scores than these. Default: %(default)s')
  cons.add_argument('-c', '--cons-thres', type=float,
    help='The threshold for calling a consensus base. More than this fraction of bases must agree '
         'in order to use the plurality vote as the consensus. E.g. give 0.5 to require a majority. '
         'Default: there is no threshold. The plurality base will be used.')
  msa = parser.add_argument_group('MSA Input')
  msa.add_argument('-b', '--barcode',
    help='The barcode of the family you want to view.')
  msa.add_argument('-o', '--order', choices=('ab', 'ba'),
    help='The order of the family you want to view.')
  msa.add_argument('-m', '--mate', choices=('0', '1', '2'),
    help='The mate of the family you want to view.')
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  log.add_argument('-s', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  log.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  log.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  format = get_format(args.format, args.input)

  seq_col, qual_col = get_columns(args.seq_column, args.qual_column, format, args.format)

  if args.format == 'msa':
    input = filter_msa(args.input, args.barcode, args.order, args.mate)
  else:
    input = args.input

  seqs, quals, seqlen = read_seqs(input, format, args.qual_format, seq_col, qual_col)

  if len(seqs) == 0 or len(quals) == 0:
    fail('Error: No sequences found.')

  masked_seqs, consensus = mask_seqs(seqs, quals, seqlen, args.qual_thres, args.cons_thres)

  print(consensus)
  for seq in masked_seqs:
    print(seq)


def get_format(format, input):
  if format:
    if format == 'msa':
      return 'tsv'
  else:
    if input is sys.stdin:
      fail('Error: Must give the --format if reading from stdin.')
    ext = os.path.splitext(input.name)[1]
    if ext == '.fq':
      return 'fastq'
    elif ext == '.fa':
      return 'fasta'
    else:
      return ext[1:]


def get_columns(seq_column, qual_column, format, raw_format):
  if raw_format == 'msa':
    seq_col = 5
    qual_col = 6
  else:
    seq_col = DEFAULT_SEQ_COLUMN
    qual_col = DEFAULT_QUAL_COLUMN
  if seq_column or qual_column:
    if format != 'tsv':
      fail('Error: --seq-column and --qual-column can only be used with tsv format.')
    if seq_column:
      seq_col = seq_column
    if qual_column:
      qual_col = qual_column
    else:
      qual_col = None
  return seq_col, qual_col


def filter_msa(input, barcode, order, mate):
  for line in input:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 3:
      return
    if barcode is not None and fields[0] != barcode:
      continue
    if order is not None and fields[1] != order:
      continue
    if mate is not None and fields[2] != mate:
      continue
    yield line


def read_seqs(infile, format, qual_format, seq_col, qual_col):
  seqlen = None
  seqs = []
  quals = []
  for read in getreads.getparser(infile, format, qual_format=qual_format,
                                 seq_col=seq_col, qual_col=qual_col):
    if not read.seq:
      continue
    if seqlen is None:
      seqlen = len(read.seq)
    if seqlen != len(read.seq):
      fail('Error: Line lengths not equal ({0} != {1})'.format(seqlen, len(read.seq)))
    if read.qual and len(read.seq) != len(read.qual):
      fail('Error: Sequence and quality scores not the same length.')
    seqs.append(read.seq)
    quals.append(read.scores)
  return seqs, quals, seqlen


def read_quals(infile, seqlen, offset):
  all_quals = []
  for line_num, line_raw in enumerate(infile):
    line = line_raw.rstrip('\r\n')
    if len(line) != seqlen:
      fail('Error: Quality scores on line {} are a different length than the input alignment '
           '({} != {}).'.format(line_num+1, len(line), seqlen))
    quals = [ord(qual) - offset for qual in line]
    all_quals.append(quals)
  return all_quals


def mask_seqs(seqs, quals, seqlen, qual_thres, cons_thres):
  mismatches = 0
  consensus = ''
  masked_seqs = [''] * len(seqs)
  for b in range(seqlen):
    votes = collections.defaultdict(int)
    for s in range(len(seqs)):
      if good_quality(s, b, quals, qual_thres):
        votes[seqs[s][b]] += 1
    max_vote = 0
    cons_base = 'N'
    for base, vote in votes.items():
      # The consensus base must be in the majority (>50%) or no consensus.
      if vote > max_vote and (cons_thres is None or vote > len(seqs)*cons_thres):
        max_vote = vote
        cons_base = base
    consensus += cons_base
    for s in range(len(seqs)):
      if not good_quality(s, b, quals, qual_thres):
        masked_seqs[s] += ' '
        continue
      if seqs[s][b] == cons_base:
        masked_seqs[s] += '.'
      else:
        mismatches += 1
        masked_seqs[s] += seqs[s][b]
  logging.info('{} mismatches'.format(mismatches))
  return masked_seqs, consensus


def good_quality(s, b, quals, qual_thres):
  if quals and quals[s]:
    q = quals[s][b]
    if 0 <= q < qual_thres:
      return False
  return True


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
