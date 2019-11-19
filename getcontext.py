#!/usr/bin/env python3
import argparse
import collections
import logging
import pathlib
import sys
from typing import Optional, Any, Iterable, Iterator, Sequence, Generator, Dict, List, Tuple
from fastagenerators import FastaLineBuffered
assert sys.version_info.major >= 3, 'Python 3 required'


DESCRIPTION = """"""


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  io = parser.add_argument_group('I/O')
  io.add_argument('ref', type=argparse.FileType('r'),
    help='Reference sequence.')
  io.add_argument('sites', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='File containing the coordinates of the sites to get sequence context for. '
      'Should be tab-delimited. Default: stdin.')
  options = parser.add_argument_group('Options')
  options.add_argument('-f', '--field', type=int, default=1,
    help='Which column in the sites file has the coordinates of the sites. Numbers are 1-based.')
  options.add_argument('-c', '--chrom-field', type=int,
    help='Which column in the sites file has the id of the chromosome of the site. '
      "Numbers are 1-based. Can omit if there's only one chromosome.")
  options.add_argument('-C', '--chrom-id',
    help='Assume all sites are in this reference sequence.')
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


def main(argv: List[str]) -> int:

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.chrom_field and args.chrom_id:
    parser.print_help()
    fail('--chrom-field and --chrom-id are mutually exclusive.')
  elif not (args.chrom_field or args.chrom_id):
    parser.print_help()
    fail('Must provide a --chrom-field or a --chrom-id.')

  sites_by_chrom = read_sites(args.sites, args.field, args.chrom_field, args.chrom_id)

  if len(sites_by_chrom) == 0:
    fail('No sites found in input.')

  # for chrom_id, sites in sites_by_chrom.items():
  #   for coord in sites:
  #     print(chrom_id, coord, sep='\t')

  for data in read_ref(args.ref, sites_by_chrom):
    print(*data, sep='\t')
    pass

  return 0


def read_sites(
    infile: Iterable[str], coord_col: int, chrom_col: int=None, chrom_id=None, sep='\t'
  ) -> Dict[str,List[int]]:
  sites_by_chrom: collections.defaultdict = collections.defaultdict(list)
  const_chrom_id = chrom_id
  chrom_id = None
  warned_value = False
  warned_index = False
  for line in infile:
    fields = line.rstrip('\r\n').split(sep)
    chrom_id, coord_str, warned_index = get_field_strs(
      fields, coord_col, chrom_col, const_chrom_id, line, warned_index
    )
    if coord_str is None:
      continue
    coord = parse_coord_int(coord_str, coord_col, warned_value)
    if coord is None:
      continue
    sites_by_chrom[chrom_id].append(coord)
  # Sort the lists of coordinates.
  for sites in sites_by_chrom.values():
    sites.sort()
  return sites_by_chrom


def get_field_strs(
    fields: List[str], coord_col: int, chrom_col: Optional[int], const_chrom_id: Optional[str],
    line: str, warned_index: bool
  ):
  coord_str = chrom_id = None
  try:
    coord_str = fields[coord_col-1]
    if const_chrom_id:
      chrom_id = const_chrom_id
    elif chrom_col is not None:
      chrom_id = fields[chrom_col-1]
  except IndexError as error:
    if not warned_index:
      logging.warning(
        f'Warning: Line(s) with not enough fields in sites file. Looking for columns '
        f'{coord_col} and {chrom_col}. Example invalid line: {line!r}'
      )
      warned_index = True
    return None, None, warned_index
  return chrom_id, coord_str, warned_index


def parse_coord_int(coord_str: str, coord_col: int, warned_value: bool):
  try:
    coord = int(coord_str)
  except ValueError:
    if not warned_value:
      logging.warning(
        f'Warning: Line(s) with invalid coordinate in sites file. Example: Found '
        f'{coord_str!r} in column {coord_col}.'
      )
      warned_value = True
    return None
  return coord


def read_ref(ref_path: pathlib.Path, sites_by_chrom: Dict[str,List[int]]):
  fasta = FastaLineBuffered(ref_path)
  sites: Iterator[int]
  site: Optional[int]
  for chrom in fasta:
    if chrom.id in sites_by_chrom and sites_by_chrom[chrom.id]:
      sites = iter(sites_by_chrom[chrom.id])
    else:
      continue
    site = next(sites)
    for coord, base in enumerate(chrom.bases(), 1):
      while coord == site:
        yield chrom.id, coord, base
        try:
          site = next(sites)
        except StopIteration:
          site = None
      if site is None:
        break
    remaining = len(list(sites))
    if site is not None:
      remaining += 1
    if remaining:
      logging.warning(f'Warning: {remaining} site(s) not found in sequence {chrom.id}')


def fail(message: Any) -> None:
  logging.critical('Error: '+str(message))
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
