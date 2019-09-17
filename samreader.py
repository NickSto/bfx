# A very simple module to hold code for quick parsing of SAM file fields
import re
try:
  import cigarlib
except ImportError:
  from bfx import cigarlib

__version__ = '0.8'
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

class Alignment(object):
  def __init__(self, fields=None, **kwargs):
    if fields is not None:
      if len(fields) < 11:
        raise FormatError('Fewer than 11 fields (only saw {})'.format(len(fields)))
      for i, field_info in enumerate(FIELDS):
        name = field_info['name']
        if fields[i] == NULL_STR:
          setattr(self, name, None)
        else:
          try:
            value = field_info['type'](fields[i])
          except ValueError:
            raise FormatError('Invalid value in {} column: {!r}'.format(name.upper(), fields[i]))
          setattr(self, name, value)
      self.tags = fields[11:]
    if kwargs:
      for field_info in FIELDS:
        name = field_info['name']
        if name in kwargs:
          setattr(self, name, kwargs[name])
      if 'tags' in kwargs:
        self.tags = tags
    self._length = None
  # Okay, here's the only non-simple bit:
  @property
  def length(self):
    if self._length is None:
      self._length = self._compute_read_length(self.seq, self.cigar)
    return self._length
  @classmethod
  def _compute_read_length(cls, seq, cigar):
    """Compute the read length from the length of the SEQ field plus any hard-clipped bases."""
    length = len(seq)
    if cigar is None:
      return length
    actions = cigarlib.split_cigar(cigar)
    for oplen, op in actions:
      if op == 'H':
        length += oplen
    return length


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
  for line_raw in filehandle:
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
    yield Alignment(fields)


def is_header(fields):
  return fields[0].startswith('@') and re.search(HEADER_REGEX, fields[0])


class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)
