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

def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('input', metavar='align.fa', type=argparse.FileType('r'), default=sys.stdin,
                      nargs='?',
    help='Aligned input sequences. Can be FASTA or just the plain sequences, one per line.')
  parser.add_argument('-f', '--format', choices=getreads.FORMATS,
    help='Format of the input. Will be inferred from the filename if not given. "lines" is the most '
         'basic format: Just the sequences, one per line.')
  parser.add_argument('-q', '--qual-thres', type=int, default=25,
    help='Quality score threshold. If quality scores are present, don\'t show bases with lower '
         'quality scores than these. Default: %(default)s')
  parser.add_argument('-F', '--qual-format', choices=QUAL_OFFSETS.keys(), default='sanger',
    help='FASTQ quality score format. Sanger scores are assumed to begin at \'{}\' ({}). '
         'Default: %(default)s.'.format(QUAL_OFFSETS['sanger'], chr(QUAL_OFFSETS['sanger'])))
  parser.add_argument('-c', '--cons-thres', type=float,
    help='The threshold for calling a consensus base. More than this fraction of bases must agree '
         'in order to use the plurality vote as the consensus. E.g. give 0.5 to require a majority. '
         'Default: there is no threshold. The plurality base will be used.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  if args.format:
    format = args.format
  else:
    if args.input is sys.stdin:
      fail('Error: Must give the --format if reading from stdin.')
    ext = os.path.splitext(args.input.name)[1]
    if ext == '.fq':
      format = 'fastq'
    elif ext == '.fa':
      format = 'fasta'
    else:
      format = ext[1:]

  seqs, quals, seqlen = read_seqs(args.input, format, args.qual_format)

  if len(seqs) == 0 or len(quals) == 0:
    fail('Error: No sequences found.')

  masked_seqs, consensus = mask_seqs(seqs, quals, seqlen, args.qual_thres, args.cons_thres)

  print(consensus)
  for seq in masked_seqs:
    print(seq)


def read_seqs(infile, format, qual_format):
  seqlen = None
  seqs = []
  quals = []
  for read in getreads.getparser(infile, format, qual_format=qual_format):
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
