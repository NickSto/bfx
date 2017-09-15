#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse
import collections
assert sys.version_info.major >= 3, 'Python 3 required'

DESCRIPTION = """Take an alignment and print a version with a consensus sequence and the conserved
bases masked."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('input', metavar='align.fa', type=argparse.FileType('r'), default=sys.stdin,
    help='Aligned input sequences. Can be FASTA or just the plain sequences, one per line.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  seqs, seqlen = read_seqs(args.input)

  masked_seqs, consensus = mask_seqs(seqs, seqlen)

  print(consensus)
  for seq in masked_seqs:
    print(seq)


def read_seqs(infile):
  seqlen = None
  seqs = []
  for line_raw in infile:
    line = line_raw.rstrip('\r\n')
    if line.startswith('>'):
      continue
    else:
      if seqlen is None:
        seqlen = len(line)
      elif seqlen != len(line):
        fail('Error: Line lengths not equal ({0} != {1})'.format(seqlen, len(line)))
      seqs.append(line)
  return seqs, seqlen


def mask_seqs(seqs, seqlen):
  mismatches = 0
  consensus = ''
  masked_seqs = [''] * len(seqs)
  for b in range(seqlen):
    votes = collections.defaultdict(int)
    for s in range(len(seqs)):
      votes[seqs[s][b]] += 1
    max_vote = 0
    cons_base = 'N'
    for base, vote in votes.items():
      if vote > max_vote:
        max_vote = vote
        cons_base = base
    consensus += cons_base
    for s in range(len(seqs)):
      if seqs[s][b] == cons_base:
        masked_seqs[s] += '.'
      else:
        mismatches += 1
        masked_seqs[s] += seqs[s][b]
  logging.info('{} mismatches'.format(mismatches))
  return masked_seqs, consensus


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
