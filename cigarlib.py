#!/usr/bin/env python3
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division
import argparse
import logging
import re
import sys

DESCRIPTION = """All the fiddly CIGAR-dependent calculations in dealing with reads."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('input', metavar='align.sam', nargs='?')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):
  import samreader

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.input:
    if args.input.endswith('.bam'):
      infile = open_bam(args.input)
    else:
      infile = open(args.sam)
  else:
    infile = sys.stdin

  for align in samreader.read(infile):
    print('pos: {}\tcigar: {}\tseq: {}'.format(align.pos, align.cigar, align.seq))
    cigar_list = split_cigar(align.cigar)
    reverse = align.flag & 16
    blocks = get_contiguous_blocks(align.pos, cigar_list, reverse, len(align.seq))
    for block in blocks:
      print('  {}'.format(block))


def open_bam(bam_path):
  import subprocess
  process = subprocess.Popen(('samtools', 'view', bam_path), stdout=subprocess.PIPE)
  if sys.version_info.major <= 2:
    return process.stdout
  else:
    return (str(line, 'utf8') for line in process.stdout)


def split_cigar(cigar):
  cigar_list = []
  bits = re.findall(r'\d+[A-Z=]', cigar)
  for bit in bits:
    op_char = bit[-1:]
    length = int(bit[:-1])
    cigar_list.append((length, op_char))
  return cigar_list


def to_ref_coord(blocks, read_coord):
  for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
    if direction == 1:
      hit = read_start <= read_coord < read_end
    elif direction == -1:
      hit = read_end < read_coord <= read_start
    if hit:
      return direction * read_coord + offset
  return None
  # logging.warn('No hit on read coordinate {}.'.format(read_coord))


#TODO: def to_read_coord(blocks, ref_pos):


def get_contiguous_blocks(ref_pos, cigar_list, reverse, read_len):
  """Return a list of blocks of aligned bases.
  Each block is a 6-tuple:
  1. The start of the block, in read coordinates.
  2. The end of the block, in read coordinates.
  3. The start of the block, in reference coordinates.
  4. The end of the block, in reference coordinates.
  5. The offset between read and reference coordinates (offset = ref - direction * read).
  6. The direction of the read (1 for forward, -1 for reverse)."""
  ref_pos_start = ref_pos
  if reverse:
    direction = -1
    read_pos = read_len
  else:
    direction = 1
    read_pos = 1
  read_pos_start = read_pos
  blocks = []
  # logging.info('Ref starting at {}, read {}.'.format(ref_pos, read_pos))
  while cigar_list:
    cigar_size, cigar_op = cigar_list.pop( 0 )
    # logging.info('Saw {}{}.'.format(cigar_size, OP_INTS_TO_CHARS[cigar_op]))
    # M alignment match (can be a sequence match or mismatch)
    # = sequence match
    # X sequence mismatch
    if cigar_op in 'M=X':
      ref_pos += cigar_size
      read_pos += cigar_size * direction
    # I insertion
    # S soft clipping (clipped sequences present in SEQ)
    elif cigar_op in 'IS':
      offset = ref_pos - direction * read_pos
      if read_pos_start != read_pos:
        blocks.append((read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction))
      read_pos += cigar_size * direction
      read_pos_start = read_pos
      ref_pos_start = ref_pos
    # D deletion from the reference
    # N skipped region from the reference
    elif cigar_op in 'DN':
      offset = ref_pos - direction * read_pos
      blocks.append((read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction))
      ref_pos += cigar_size
      read_pos_start = read_pos
      ref_pos_start = ref_pos
    # H hard clipping (clipped sequences NOT present in SEQ)
    # P padding (silent deletion from padded reference)
    elif cigar_op in 'HP':
      pass
    else:
      pass #logging.warn('unknown cigar_op {} {}'.format(cigar_op, cigar_size))
    # logging.info('Ref now {}, read {}.'.format(ref_pos, read_pos))
  offset = ref_pos - direction * read_pos
  if read_pos_start != read_pos:
    blocks.append((read_pos_start, read_pos, ref_pos_start, ref_pos, offset, direction))
  return blocks


def indel_at(position, insertions, deletions, check_insertions=True, check_deletions=True):
  """Does the read contain an indel at the given position?
  Return True if the read contains an insertion at the given position
  (position must be the base before the insertion event) or if the read
  contains a deletion where the base at position is deleted. Return False
  otherwise."""
  if check_insertions:
    for insertion in insertions:
      if insertion[0] == position:
        return True
  if check_deletions:
    for deletion in deletions:
      if deletion[0] < position < deletion[0] + deletion[1] + 1:
        return True
  return False


def get_indels(blocks, reverse):
  """Return a data structure containing all indels in the read.
  Returns the tuple (insertions, deletions)
  insertions = [(pos1,ins1), (pos2,ins2)]
  posN = start position (preceding base, VCF-style)
  insN = length of inserted sequence (not including preceding base)
  deletions = [(pos1,del1), (pos2,del2)]
  posN = start position (preceding base, VCF-style)
  delN = length of deleted sequence (not including preceding base)
  Note: This does not count any "I" or "D" CIGAR operations at the start or end of a read.
  It also counts "N" as a deletion."""
  #TODO: Include the cigar operation as a field in the block so this can avoid counting "N"s
  #      as deletions.
  insertions = []
  deletions = []
  last_read_end = None
  last_ref_end = None
  if reverse:
    for read_end, read_start, ref_end, ref_start, offset, direction in reversed(blocks):
      if last_read_end is not None:
        if read_start == last_read_end:
          del_len = last_ref_end-ref_start
          del_start = last_ref_end-del_len-1
          deletions.append((del_start, del_len))
          logging.info('DEL: last_ref_end-del_len-1 = {}-{}-1 = {}; last_ref_end-ref_start = {}-{} = {}'
                       .format(last_ref_end, del_len, del_start, last_ref_end, ref_start, del_len))
        else:
          ins_start = last_ref_end-1
          ins_len = read_start-last_read_end
          insertions.append((ins_start, ins_len))
      last_read_end = read_end
      last_ref_end = ref_end
  else:
    for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
      if last_read_end is not None:
        if read_start == last_read_end:
          del_start = last_ref_end-1
          del_len = ref_start-last_ref_end
          deletions.append((del_start, del_len))
          logging.info('DEL: last_ref_end-1 = {}-1 = {}; ref_start-last_ref_end = {}-{} = {}'
                       .format(last_ref_end, del_start, ref_start, last_ref_end, del_len))
        else:
          ins_start = last_ref_end-1
          ins_len = read_start-last_read_end
          insertions.append((ins_start, ins_len))
          logging.info('INS: last_ref_end-1 = {}-1 = {}; read_start-last_read_end = {}-{} = {}'
                       .format(last_ref_end, ins_start, read_start, last_read_end, ins_len))
      last_read_end = read_end
      last_ref_end = ref_end
  return insertions, deletions


def get_end_position(blocks):
  """Get the position of the end of the aligned portion of the sequence.
  This will always be the "right-most" position, regardless of orientation of the read.
  So for a read with its 5' end at 100 and its 3' end at 200, this will give 200.
  And for a reverse-oriented read with its 3' end at 100 and its 5' end at 200, this will
  still give 200.
  Note: hard- and soft-clipped bases are not counted."""
  for read_start, read_end, ref_start, ref_end, offset, direction in blocks:
    max_position = max(ref_start, ref_end)
  return max_position


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  sys.exit(main(sys.argv))
