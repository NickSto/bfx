# A very simple module to hold code for quick parsing of SAM file fields
import re
__version__ = '0.7'
HEADER_REGEX = r'^@[A-Za-z][A-Za-z]$'
FIELD_NAMES = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext',
  'pnext', 'tlen', 'seq', 'qual')

class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)

def read(filehandle, header=False):
  """A simple SAM parser.
  This returns a generator, so it only reads a line at a time.
  Initialize with an open file (or sys.stdin).
  Yields a dict for each read, mapping the SAM field names (lowercase) to
  their values. The tags are under the key 'tags', in a list.
  If "header" argument is True, it will first yield the header lines, one at
  a time, as simple lists of the tab-delimited field values.
  Usage:
  for fields in SamReader(open('/home/user/alignment.sam', 'rU')):
    print "QNAME: "+fields['qname']+", POS: "+fields['pos']
    print len(fields['tags']), "tags present."
  """
  for line_raw in filehandle:
    line = line_raw.strip()
    # skip empty lines (whitespace-only)
    if len(line) == 0:
      continue
    fields = line.split('\t')
    # skip header lines by default
    if fields[0].startswith('@') and re.search(HEADER_REGEX, fields[0]):
      if header:
        yield fields
      else:
        continue
    if len(fields) < 11:
      raise FormatError('Fewer than 11 fields on line:\n'+line_raw)
    yield _dictify(fields, FIELD_NAMES)

def _dictify(fields, field_names):
  data = dict(zip(field_names, fields[:11]))
  data['tags'] = fields[11:]
  return data
