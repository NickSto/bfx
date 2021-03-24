#!/usr/bin/env python3
import argparse
import collections
import logging
import pathlib
import sys
import samreader

DESCRIPTION = """Get the mapped read depth at each reference position from an alignment.
WARNING: Currently cannot handle alignments to reference sequences with multiple chromosomes."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('alignment', type=pathlib.Path, nargs='?',
    help='The input alignment. Give the path to a SAM or BAM file, or omit to read SAM from stdin. '
      'The format will be inferred from the filename.')
  options.add_argument('-f', '--format', choices=('sam', 'bam'),
    help='The format of the input alignment file. This overrides any inferred format from the '
      'filename.')
  options.add_argument('-x', '--exclude', type=csv_ints, default=(),
    help='Exclude alignments with these SAM flags set. Give a comma-separated list of integers.')
  options.add_argument('-i', '--include', type=csv_ints, default=(),
    help="Exclude alignments that DON'T have these SAM flags set. Give a comma-separated list of "
      'integers.')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  align_lines = open_alignment(args.alignment, args.format)

  alignments = filter_alignments(samreader.read(align_lines), args.exclude, args.include)

  depths = get_depths(alignments)

  for coord, depth in sorted(depths.items(), key=lambda item: item[0]):
    print(coord, depth, sep='\t')


def open_alignment(align_path, format_arg):
  if align_path:
    if format_arg:
      format_ = format_arg
    else:
      format_ = get_format(align_path)
      if format_ is None:
        fail(f'Format could not be inferred from filename {align_path.name!r}')
    if format_ == 'sam':
      return align_path.open('r')
    elif format_ == 'bam':
      return samreader.open_bam(align_path)
  else:
    return sys.stdin


def get_format(path):
  if path.suffix.lower() == '.sam':
    return 'sam'
  elif path.suffix.lower() == '.bam':
    return 'bam'
  else:
    return None


def filter_alignments(alignments, exclude_flags, include_flags):
  removed = 0
  for i, align in enumerate(alignments,1):
    if any([align.has_flag(flag) for flag in exclude_flags]):
      removed += 1
      continue
    if not all([align.has_flag(flag) for flag in include_flags]):
      removed += 1
      continue
    yield align
  logging.info(f'Info: Filtered out {removed} poor quality alignments out of {i} total.')


def get_depths(alignments):
  depths = collections.Counter()
  for i, align in enumerate(alignments,1):
    for ref_coord in get_ref_coords(align):
      depths[ref_coord] += 1
  logging.info(f'Info: Found {len(depths)} reference positions in {i} alignments.')
  return depths


def get_ref_coords(align):
  for read_coord in range(1,len(align.seq)+1):
    ref_coord = align.to_ref_coord(read_coord)
    if ref_coord is not None:
      yield ref_coord


def csv_ints(csv_str):
  return [int(val) for val in csv_str.split(',')]


def fail(message):
  logging.critical(f'Error: {message}')
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
