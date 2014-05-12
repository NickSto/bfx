# A very simple module to hold code for quick parsing of SAM file fields
import re
__version__ = '0.1'
HEADER_REGEX = r'^@[A-Za-z][A-Za-z]$'

class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


class SamReader(object):
  """A simple SAM parser.
  This returns a generator, so it only reads a line at a time.
  Initialize with an open file (or sys.stdin).
  Yields a list of the SAM fields on each iteration.
  Usage:
  for fields in SamReader(open('/home/user/alignment.sam', 'rU')):
    print "QNAME: "+fields[0]+", POS: "+fields[3]
    print str(len(fields) - 11)+" tags present."
  """

  def __init__(self, filehandle):
    self.filehandle = filehandle

  def __iter__(self):
    return self.new()

  def new(self):
    for line_raw in self.filehandle:
      line = line_raw.strip()
      # skip empty lines (whitespace-only)
      if len(line) == 0:
        continue
      fields = line.split('\t')
      # skip header lines
      if fields[0].startswith('@') and re.search(HEADER_REGEX, fields[0]):
        continue
      if len(fields) < 11:
        raise FormatError('Fewer than 11 fields on line:\n'+line_raw)
      yield fields
