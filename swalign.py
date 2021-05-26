#!/usr/bin/env python3
import ctypes
import errno
import os
import sys

# Locate the library file.
LIBFILE = 'libswalign.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

swalign = ctypes.cdll.LoadLibrary(library_path)

REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')


# C struct for ctypes
class SeqPairC(ctypes.Structure):
  _fields_ = [
    ('a', ctypes.c_char_p),
    ('alen', ctypes.c_uint),
    ('b', ctypes.c_char_p),
    ('blen', ctypes.c_uint),
  ]


# C struct for ctypes
class AlignC(ctypes.Structure):
  _fields_ = [
    ('seqs', ctypes.POINTER(SeqPairC)),
    ('start_a', ctypes.c_int),
    ('start_b', ctypes.c_int),
    ('end_a', ctypes.c_int),
    ('end_b', ctypes.c_int),
    ('matches', ctypes.c_int),
    ('score', ctypes.c_double),
  ]


# The Python version
class Align(object):
  def __init__(self, align_c):
    self.target = str(align_c.seqs.contents.a, 'utf8')
    self.query = str(align_c.seqs.contents.b, 'utf8')
    # Where the first base of the target aligns on the query, in query coordinates (or 1, if <= 0).
    self.start_target = align_c.start_a
    # Where the first base of the query aligns on the target, in target coordinates (or 1, if <= 0).
    self.start_query = align_c.start_b
    # Where the last base of the target aligns on the query, in query coordinates.
    self.end_target = align_c.end_a
    # Where the last base of the query aligns on the target, in target coordinates.
    self.end_query = align_c.end_b
    self.matches = align_c.matches
    self.score = align_c.score

  # Provide this common function.
  def __str__(self):
    """Print a human-readable representation of the alignment."""
    start_query = str(self.start_query)
    start_target = str(self.start_target)
    start_width = max(len(start_query), len(start_target))
    line_format = '{:'+str(start_width)+'} {} {}'
    output = line_format.format(start_target, self.target, self.end_target) + '\n'
    output += ' '*(start_width+1) + format_matches(self.target, self.query) + '\n'
    output += line_format.format(start_query, self.query, self.end_query)
    return output


def format_matches(target, query):
  matches_str = ''
  for base1, base2 in zip(target, query):
    if base1 == base2 and base1 != '-':
      matches_str += '|'
    else:
      matches_str += ' '
  return matches_str


# Initialize functions (define types).
swalign.smith_waterman.restype = ctypes.POINTER(AlignC)
swalign.revcomp.restype = ctypes.c_char_p


def smith_waterman(target_raw, query_raw, local=True, debug=False):
  target_bytes = bytes(target_raw, 'utf8')
  query_bytes = bytes(query_raw, 'utf8')
  seq_pair = SeqPairC(target_bytes, len(target_raw), query_bytes, len(query_raw))
  align_c = swalign.smith_waterman(ctypes.pointer(seq_pair), bool(local), bool(debug)).contents
  return Align(align_c)


def smith_waterman_duplex(target, query):
  """Smith-Waterman align query to target in both orientations and return the best.
  Convenience function that calls smith_waterman() twice, and returns the
  alignment with the highest score."""
  align = smith_waterman(target, query)
  query_rc = revcomp(query)
  align_rc = smith_waterman(target, query_rc)
  if align_rc.score > align.score:
    return align_rc
  else:
    return align


def revcomp(seq):
  """Return the reverse complement of the input sequence.
  Leaves the input string unaltered."""
  return seq.translate(REVCOMP_TABLE)[::-1]


########## Command-line interface ##########

import argparse
import logging

DESCRIPTION = """Align two sequences with Smith-Waterman."""

def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('seq1')
  parser.add_argument('seq2')
  parser.add_argument('-l', '--local', action='store_true')
  parser.add_argument('-L', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  result = smith_waterman(args.seq1, args.seq2, args.local, args.volume == logging.DEBUG)

  print(result)


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  sys.exit(main(sys.argv))
