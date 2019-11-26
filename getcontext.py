#!/usr/bin/env python3
import argparse
import collections
import logging
import pathlib
import sys
from typing import Optional, Any, Iterable, Iterator, Sequence, Generator, Dict, List, Tuple
from fastagenerators import FastaLineBuffered
from utillib import simplewrap
assert sys.version_info.major >= 3, 'Python 3 required'


DESCRIPTION = simplewrap.wrap(
"""Get info on the sequence context surrounding given sites.
Output is tab-delimited, one line per site.
Columns:
1. Sequence (chromosome) name.
2. Coordinate in the sequence.
3. Coordinate in the context excerpt (0-based).
4. Base of the site in the context excerpt.
   - Semi-redundant: excerpt[coord] == base ($5[$3] == $4)
5. Context excerpt.
   - May be less than --window bp long.
     - I.e. if the site is within window/2 of the edge of sequence.
6. G/C content of the excerpt, as a percentage.
   - G, C, and S are counted as G/C.
   - A, T, and W are counted as A/T.
   - Case-insensitive.
   - All other bases are ignored.
   - G/C content = 100*sum(G/C)/(sum(G/C)+sum(A/T)).
   - Null content is given as '.'"""
)


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION,
    formatter_class=argparse.RawDescriptionHelpFormatter)
  io = parser.add_argument_group('I/O')
  io.add_argument('ref', type=argparse.FileType('r'),
    help='Reference sequence.')
  io.add_argument('sites', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='File containing the coordinates of the sites to get sequence context for. '
      'Should be tab-delimited. Default: stdin.')
  io.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
    help='Write output to this path instead of stdout.')
  options = parser.add_argument_group('Options')
  options.add_argument('-f', '--field', type=int, default=1,
    help='Which column in the sites file has the coordinates of the sites. Numbers are 1-based.')
  options.add_argument('-c', '--chrom-field', type=int,
    help='Which column in the sites file has the id of the chromosome of the site. '
      "Numbers are 1-based. Can omit if there's only one chromosome.")
  options.add_argument('-C', '--chrom-id',
    help='Assume all sites are in this reference sequence.')
  options.add_argument('-w', '--window', type=int, default=20,
    help='Width of the sequence context window, in bp. Will be centered on the site.')
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

  for chrom, coord, i, context in get_context(args.ref, sites_by_chrom, args.window):
    gc = get_gc(context, null='.', decimals=1)
    print(chrom, coord, i, context[i], context, gc, sep='\t', file=args.output)

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


def get_gc(seq: str, null: Any=None, decimals: int=None) -> Optional[float]:
  gc = 0
  at = 0
  for base in seq.upper():
    if base in ('G', 'C', 'S'):
      gc += 1
    elif base in ('A', 'T', 'W'):
      at += 1
  total = at+gc
  if total:
    percent = 100*gc/total
    if decimals is None:
      return percent
    else:
      return round(percent, decimals)
  else:
    return null


def get_context(
    ref_path: pathlib.Path, sites_by_chrom: Dict[str,List[int]], window: int
  ) -> Generator[Tuple[str,int,int,str],None,None]:
  fasta = FastaLineBuffered(ref_path)
  sites: Iterator[int]
  site: Optional[int]
  for chrom in fasta:
    if chrom.id in sites_by_chrom and sites_by_chrom[chrom.id]:
      sites = iter(sites_by_chrom[chrom.id])
    else:
      continue
    site = get_next_site(sites)
    context = Context(window=window)
    for coord, base in enumerate(chrom.bases(), 1):
      context.push(base)
      while site is not None and site == context.middle:
        yield chrom.id, site, context.middle_index, str(context)
        site = get_next_site(sites)
      if site is None:
        break
    while site is not None:
      context.shift()
      while site is not None and site == context.middle:
        yield chrom.id, site, context.middle_index, str(context)
        site = get_next_site(sites)
    remaining = len(list(sites))
    if site is not None:
      remaining += 1
    if remaining:
      logging.warning(f'Warning: {remaining} site(s) not found in sequence {chrom.id}')


def get_next_site(sites):
  try:
    site = next(sites)
  except StopIteration:
    site = None
  return site


class Context:
  """An object representing a sequence context using a sliding window.
  This tracks both the window and the sequence contained within it. The two are
  separate because the window can hang off the left or right edge of the sequence,
  with the left edge in negative coordinates, for example.
  The distance between the left and right edges of the window is always constant,
  i.e. `context.right - context.left + 1 == context.window` always.
  But the sequence in the `Context` can be smaller than `context.window`, for
  example if the window is past the left or right edges of the sequence.
  Note:
  All coordinates are 1-based, where each base has a coordinate, starting with 1.
  Attributes:
  `window`: The size of the window (integer).
  `left` and `right`: Coordinates of the window edges. `left` can be <= 0, and `right`
    can be greater than the sequence length.
  `middle`: Coordinate of the center of the window. If the window size is even, this
    will be the coordinate to the left of the center.
  `left_base`, `right_base`, `middle_base`: The bases at each of these coordinates.
    If the coordinate is outside the sequence bounds, the base will be `None`.
  `seq_left`, `seq_right`: The coordinates of the first and last bases contained in
    the window. These may be different than `left` and `right` if the window is past
    the edge of the sequence.
  `middle_index`: The index of the middle base in the string representation of this
    `Context`. `str(context)[context.middle_index]` should be equivalent to
    `context.middle_base`."""

  def __init__(self, seq: str=None, window=10):
    #TODO: If performance is an issue, try a list, which may be faster for small windows.
    self._deque: collections.deque = collections.deque()
    self._window = window
    self.left = 1-self.window
    self.right = 0
    self.seq_left = 1
    self.seq_right = 0
    if seq:
      for base in seq:
        self.push(base)

  def push(self, base: str) -> Optional[str]:
    """Add a new base to the right side of the context and slide the window to the right."""
    self._deque.append(base)
    self.left += 1
    self.right += 1
    self.seq_right += 1
    return self._fix_left_end()

  def shift(self):
    """Move the window right without adding more bases.
    Raises an `IndexError` once there's no more bases left in the context."""
    self.left += 1
    self.right += 1
    return self._fix_left_end()

  def _fix_left_end(self) -> Optional[str]:
    popped_base = None
    while (
        self.seq_left < self.left or
        self.seq_left < self.seq_right-self.window+1 or
        len(self) > self.window
    ):
      popped_base = self._deque.popleft()
      self.seq_left += 1
    return popped_base

  @property
  def window(self) -> int:
    return self._window

  @property
  def middle(self) -> int:
    """Return the coordinate of the "middle" base in the window.
    If the context is smaller than the window, this can return a coordinate less than
    `self.left` (even a negative coord)."""
    return self.right-(self.window//2)

  @property
  def middle_base(self) -> Optional[str]:
    middle = self.middle
    try:
      return self[middle]
    except IndexError:
      return None

  @property
  def middle_index(self):
    return self.middle - self.seq_left

  @property
  def left_base(self):
    try:
      return self[self.seq_left]
    except IndexError:
      return None

  @property
  def right_base(self):
    try:
      return self[self.seq_right]
    except IndexError:
      return None

  def __len__(self) -> int:
    return len(self._deque)

  def __iter__(self) -> Iterator[str]:
    return iter(self._deque)

  def __getitem__(self, coord: int) -> str:
    """`coord` should be in the range `self.seq_left <= coord <= self.seq_right`"""
    index = coord - self.seq_left
    base = None
    if index >= 0:
      try:
        base = self._deque[index]
      except IndexError:
        pass
    if base is None:
      raise IndexError(
        f'Coordinate {coord} out of range ({type(self).__name__} object covers bases '
        f'{self.seq_left} to {self.seq_right}).'
      )
    else:
      return base

  def __str__(self):
    return ''.join(self._deque)

  def __repr__(self):
    class_name = type(self).__name__
    return f"<{class_name} object {str(self)!r}>"


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
