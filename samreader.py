#!/usr/bin/env python3
# A not-so-simple module to hold code for quick parsing of SAM file fields
import argparse
import logging
import re
import subprocess
import sys
try:
  import cigarlib
except ImportError:
  from bfx import cigarlib

__version__ = '0.10'
NULL_STR = '*'
HEADER_REGEX = r'^@[A-Za-z][A-Za-z]$'
FIELDS = (
  {'name':'qname', 'type':str},
  {'name':'flag', 'type':int},
  {'name':'rname', 'type':str},
  {'name':'pos', 'type':int},
  {'name':'mapq', 'type':int},
  {'name':'cigar', 'type':str},
  {'name':'rnext', 'type':str},
  {'name':'pnext', 'type':int},
  {'name':'tlen', 'type':int},
  {'name':'seq', 'type':str},
  {'name':'qual', 'type':str},
)
FLAGS = {
  'paired': 1,
  'proper': 2,
  'unmapped': 4,
  'mapped':  -4,
  'mate_unmapped': 8,
  'reverse':  16,
  'reversed': 16,
  'forward': -16,
  'mate_reverse': 32,
  'first': 64,
  'second': 128,
  'secondary': 256,
  'primary':  -256,
  'lowqual': 512,
  'duplicate': 1024,
  'supplemental': 2048,
}
FLAG_DESCRIPTIONS = {
  1: 'Read paired',
  2: 'Read mapped in proper pair',
  4: 'Read unmapped',
  8: 'Mate unmapped',
  16: 'Read reverse strand',
  32: 'Mate reverse strand',
  64: 'First in pair',
  128: 'Second in pair',
  256: 'Not primary alignment',
  512: 'Read fails platform/vendor quality checks',
  1024: 'Read is pcr or optical duplicate',
  2048: 'Supplementary alignment',
}

class Alignment:
  def __init__(self, fields=None, **kwargs):
    self._cache = {}
    self._raw_tags = None
    if fields is not None:
      if len(fields) < 11:
        raise FormatError(f'Fewer than 11 fields (only saw {len(fields)})')
      for i, field_info in enumerate(FIELDS):
        name = field_info['name']
        if fields[i] == NULL_STR:
          setattr(self, name, None)
        else:
          try:
            value = field_info['type'](fields[i])
          except ValueError:
            raise FormatError(f'Invalid value in {name.upper()} column: {fields[i]!r}')
          setattr(self, name, value)
      self._raw_tags = fields[11:]
    if kwargs:
      for field_info in FIELDS:
        name = field_info['name']
        if name in kwargs:
          setattr(self, name, kwargs[name])
      if 'tags' in kwargs:
        self.tags = kwargs['tags']
      self.line_num = kwargs.get('line_num')
  @property
  def mate(self):
    if self.has_flag(64):
      return 1
    elif self.has_flag(128):
      return 2
  # Present an attribute for each flag bit (and some of their opposites).
  def __getattr__(self, attr):
    try:
      flag_int = FLAGS[attr]
      if flag_int >= 0:
        return self.has_flag(flag_int)
      else:
        return not self.has_flag(-flag_int)
    except KeyError:
      raise AttributeError(f'{type(self).__name__!r} object has no attribute named {attr!r}')
  def has_flag(self, bit):
    return bool(self.flag & bit)
  def explain_flag(self, no_print=False):
    outlines = []
    for flag_int, description in FLAG_DESCRIPTIONS.items():
      if self.has_flag(flag_int):
        check = 'X'
      else:
        check = ' '
      outlines.append(f'[{check}] {flag_int:4d} {description}')
    output = '\n'.join(outlines)
    if no_print:
      return output
    else:
      print(output)
  # Read length.
  @property
  def length(self):
    if 'length' not in self._cache:
      if self.cigar is None:
        #TODO: Maybe could just return len(self.seq). If it's not mapped, then there presumably
        #      wouldn't be any clipping. But would need to double-check that assumption.
        raise FormatError(
          f'No CIGAR string present. Cannot determine read length.', line_num=self.line_num
        )
      self._cache['length'] = self._compute_read_length(self.seq, self._cigar_list)
    return self._cache['length']
  @classmethod
  def _compute_read_length(cls, seq, cigar_list):
    """Compute the read length from the length of the SEQ field plus any hard-clipped bases."""
    length = len(seq)
    for oplen, op in cigar_list:
      if op == 'H':
        length += oplen
    return length
  # Tags.
  @property
  def tags(self):
    self._ensure_tags_are_parsed()
    return self._cache['tags']
  @tags.setter
  def tags(self, value):
    if not isinstance(value, dict):
      raise ValueError(f'Tags must be a dict. Received a {type(value)} instead.')
    self._cache['tags'] = value
  @property
  def tag_types(self):
    self._ensure_tags_are_parsed()
    return self._cache['tag_types']
  @tag_types.setter
  def tag_types(self, value):
    if not isinstance(value, dict):
      raise ValueError(f'Tag types must be a dict. Received a {type(value)} instead.')
    self._cache['tag_types'] = value
  def _ensure_tags_are_parsed(self):
    # Parse tags lazily.
    if 'tags' not in self._cache or 'tag_types' not in self._cache:
      try:
        self._cache['tags'], self._cache['tag_types'] = self._parse_tags(self._raw_tags)
      except FormatError as error:
        error.line_num = self.line_num
        raise error
  @classmethod
  def _parse_tags(cls, tags_list):
    tags = {}
    types = {}
    if tags_list is None:
      return tags, types
    for tag_str in tags_list:
      tag_name, tag_type, value = cls._parse_tag(tag_str)
      if tag_name in tags:
        raise FormatError(f'Multiple {tag_name!r} tags seen in one line.')
      tags[tag_name] = value
      types[tag_name] = tag_type
    return tags, types
  @classmethod
  def _parse_tag(cls, tag_str):
    fields = tag_str.split(':')
    if len(fields) < 3:
      raise FormatError(f'Invalid tag. Saw {len(fields)} fields instead of 3: {tag_str!r}')
    tag_name = fields[0]
    tag_type = fields[1]
    value_str = ':'.join(fields[2:])
    if not (
      len(tag_name) == 2 and (tag_name[0].isupper() or tag_name[0].islower()) and
      (tag_name[1].isupper() or tag_name[1].islower() or tag_name[1].isnumeric())
    ):
      raise FormatError(f'Invalid tag name {tag_name!r} in {tag_str!r}')
    value = cls._cast_tag(tag_name, tag_type, value_str)
    return tag_name, tag_type, value
  @classmethod
  def _cast_tag(cls, tag_name, tag_type, value_str):
    if tag_type == 'A':
      if len(value_str) != 1:
        raise make_tag_error(tag_name, tag_type, value_str, 'value is not a single character')
      if value_str in ' \t\r\n':
        raise make_tag_error(tag_name, tag_type, value_str, 'value is whitespace')
      value = value_str
    elif tag_type == 'i':
      try:
        value = int(value_str)
      except ValueError:
        raise make_tag_error(tag_name, tag_type, value_str, 'value is not an integer')
    elif tag_type == 'f':
      try:
        value = float(value_str)
      except ValueError:
        raise make_tag_error(tag_name, tag_type, value_str, 'value is not a float')
    elif tag_type == 'Z':
      if len(value_str) < 1:
        raise make_tag_error(tag_name, tag_type, value_str, 'value is empty')
      value = value_str
    elif tag_type == 'H':
      try:
        value = bytes.fromhex(value_str)
      except ValueError:
        raise make_tag_error(tag_name, tag_type, value_str, 'value is not valid hex')
    elif tag_type == 'B':
      if not (len(value_str) >= 1 and value_str[0] in '[cCsSiIf]'):
        raise make_tag_error(tag_name, tag_type, value_str, 'value does not contain a valid type')
      #TODO: Properly parse array.
      value = value_str
    elif tag_type == 'S':
      value = value_str
    else:
      tag_str = f'{tag_name}:{tag_type}:{value_str}'
      raise FormatError(f'Invalid tag. Unrecognized type {tag_type!r} in {tag_str!r}')
    return value
  # CIGAR stuff. Basically just providing direct (and cached) access to cigarlib.
  @property
  def _cigar_list(self):
    if 'cigar_list' not in self._cache:
      self._cache['cigar_list'] = cigarlib.split_cigar(self.cigar)
    return self._cache['cigar_list']
  @property
  def _contiguous_blocks(self):
    if 'contiguous_blocks' not in self._cache:
      self._cache['contiguous_blocks'] = cigarlib.get_contiguous_blocks(
        self.pos, self._cigar_list, self.reverse, len(self.seq)
      )
    return self._cache['contiguous_blocks']
  @property
  def _insertions(self):
    self._ensure_indels_are_parsed()
    return self._cache['insertions']
  @property
  def _deletions(self):
    self._ensure_indels_are_parsed()
    return self._cache['deletions']
  def _ensure_indels_are_parsed(self):
    if 'insertions' not in self._cache or 'deletions' not in self._cache:
      indels = cigarlib.get_indels(self._contiguous_blocks, self.reverse)
      self._cache['insertions'], self._cache['deletions'] = indels
  def indel_at(self, position, check_insertions=True, check_deletions=True):
    return cigarlib.indel_at(
      position, self._insertions, self._deletions, check_insertions=True, check_deletions=True
    )
  def to_ref_coord(self, read_coord):
    if self.mapped:
      return cigarlib.to_ref_coord(self._contiguous_blocks, read_coord)
    else:
      return None
  def get_end_position(self):
    return cigarlib.get_end_position(self._contiguous_blocks)
  def __repr__(self):
    return f'<{type(self).__name__} {self.qname}/{self.mate}>'


def make_tag_error(tag_name, tag_type, value_str, message):
  return FormatError(f'Invalid {tag_name} tag. Type is {tag_type!r} but {message}: {value_str!r}')


def read_bam(bam_path, header=False):
  sam = open_bam(bam_path)
  return read(sam, header=header)


def open_bam(bam_path):
  command = ('samtools', 'view', '-h', bam_path)
  process = subprocess.Popen(command, encoding='utf8', stdout=subprocess.PIPE)
  return process.stdout


def read(filehandle, header=False):
  """A simple SAM parser.
  This returns a generator, so it only reads a line at a time.
  Initialize with an open file (or sys.stdin).
  Yields an `Alignment` object for each read, where the values in the SAM columns are stored
  as the attributes of the object. The attribute are named after the field names (lowercase) from
  the spec.
  The tags are in the attribute `tags`, in a list.
  If "header" argument is True, it will first yield the header lines, one at
  a time, as simple lists of the tab-delimited field values.
  Usage:
  for fields in read(open('alignment.sam', 'rU')):
    print("QNAME: "+fields.qname.+", POS: "+str(fields.pos))
    print(len(fields.tags), "tags present.")
  """
  for line_num, line_raw in enumerate(filehandle, 1):
    line = line_raw.strip()
    # skip empty lines (whitespace-only)
    if len(line) == 0:
      continue
    fields = line.split('\t')
    # skip header lines by default
    if is_header(fields):
      if header:
        yield fields
      continue
    try:
      yield Alignment(fields, line_num=line_num)
    except FormatError as error:
      error.line_num = line_num
      raise error


def is_header(fields):
  return fields[0].startswith('@') and re.search(HEADER_REGEX, fields[0])


class FormatError(Exception):
  def __init__(self, message=None, line_num=None):
    self.message = message
    if message:
      Exception.__init__(self, message)
    self.line_num = line_num
  def __str__(self):
    if self.line_num is not None:
      return f'on line {self.line_num}: {self.message}'
    else:
      return self.message


#################### Command Line Interface ####################

DESCRIPTION = """Tools for working with SAM files."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('command', choices=('validate',),
    help='The action to take.')
  options.add_argument('sam', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
    help='The input SAM file. Omit to read from stdin.')
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

  if args.command == 'validate':
    validate(args.sam)


def validate(sam_file):
  total = 0
  for alignment in read(sam_file):
    len(alignment.tags)
    total += 1
  print(f'Saw {total} valid alignments.')


def fail(message):
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
