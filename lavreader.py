#!/usr/bin/env python
# requires Python 2.7
__version__ = '9d848d0'


class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


class LavReader(object):
  """
  Parse an LAV file and provide an API for querying its fields.
  By default, the data will be represented exactly as it appears in the file
  (but with numbers as int types). To convert coordinates to a more intuitive
  system, call the convert() method.
    Data structures:
  LavReader.hits        A list of LavHits
  LavHit.parent         The LavReader containing the LavHit
  LavHit.subject        A dict containing data on the hit's subject sequence
  LavHit.query          A dict containing data on the hit's query sequence
      keys for subject and query:
      filename, id, name, seqnum, revcomp, begin, end 
  LavHit.alignments     A list of LavAlignments
  LavAlignment.parent   The LavHit containing the LavAlignment
  LavAlignment.score    An int: the score of the hit
  LavAlignment.subject  A dict containing the subject start and end coordinate
  LavAlignment.query    A dict containing the query start and end coordinate
      keys for subject and query:
      begin, end
  LavAlignment.blocks   A list of LavBlocks, one for each gap-free block
  LavBlock.parent       The LavAlignment containing the LavBlock
  LavBlock.identity     An int: the percent identity of the block
  LavBlock.subject      A dict containing the subject start and end coordinate
  LavBlock.query        A dict containing the query start and end coordinate
      keys for subject and query:
      begin, end
    __len__ implementations:
  len(LavHit)           = number of alignments in LavHit.alignments
  len(LavAlignment)     = number of blocks in LavAlignment.blocks
  len(LavBlock)         = number of bases aligned in the block
    __iter__ implementations:
  iter(LavReader)       = iter(LavReader.hits)
  iter(LavHit)          = iter(LavHit.alignments)
  iter(LavAlignment)    = iter(LavAlignment.blocks)

  Format assumptions:
  **This was written only to parse the output of LASTZ version 1.02.00.**
  Stanza starts and ends are on their own lines.
  - E.g. there will be nothing (except whitespace) before or after "a {" on
    the line in which it appears. The same goes for "}".
  Stanza labels are single alphabetic characters ("Census" stanzas are ignored).
  Sequence file names do not contain whitespace.
  h stanzas are present.
  s stanzas:
  - rev_comp_flag's and sequence_number's are given.
  - only query sequences can be reverse complemented.
  """
  def __iter__(self):
    return iter(self.hits)

  def __len__(self):
    return len(self.hits)

  def __init__(self, filepath):
    self.hits = []
    self.converted = False
    
    stanza = ''
    stanza_done = True
    stanza_line = 0
    line_num = 0
    current_hit = LavHit(self)
    current_alignment = LavAlignment(current_hit)
    with open(filepath, 'rU') as filehandle:
      for raw_line in filehandle:
        line = raw_line.strip()
        line_num+=1
        if not line:
          continue
        # at start of another stanza?
        if stanza_done:
          stanza = self._stanza_start(line)
          if stanza:
            stanza_done = False
            stanza_line = 0
          continue
        # at end of a stanza?
        if line == '}':
          # if end of an alignment block, add it to the current hit
          if stanza == 'a' and len(current_alignment) > 0:
            current_hit.alignments.append(current_alignment)
            current_alignment = LavAlignment(current_hit)
          stanza_done = True
          stanza = ''
          continue
        # parse the stanzas
        stanza_line+=1
        if stanza == 's':
          # add previous hit to the list, start a new one
          if len(current_hit) > 0:
            self.hits.append(current_hit)
            current_hit = LavHit(self)
            current_alignment = LavAlignment(current_hit)
          # file names, input sequence info (revcomp, etc)
          try:
            current_hit = self._parse_s(line, current_hit, stanza_line)
          except FormatError as fe:
            fe.args = (fe.args[0]+' (line '+str(line_num)+')',)
            raise fe
        elif stanza == 'h':
          # sequence names
          current_hit = self._parse_h(line, current_hit, stanza_line)
        elif stanza == 'a':
          try:
            current_alignment = self._parse_a(line, current_alignment, stanza_line)
          except FormatError as fe:
            fe.args = (fe.args[0]+' (line '+str(line_num)+')',)
            raise fe
          except ValueError:
            raise FormatError('Invalid LAV: Non-integer encountered in "a" '
              +'stanza on line '+str(line_num)+'.')
    if len(current_alignment) > 0:
      current_hit.alignments.append(current_alignment)
    if len(current_hit) > 0:
      self.hits.append(current_hit)


  def convert(self):
    """Convert all coordinates to forward strand, 1-based, from the start.
    Also reverses the order of alignments and blocks in hits with reverse-
    complemented query sequences.
    There cannot be any reverse-complemented subject sequences, or this will
    raise a FormatError."""
    if self.converted:
      return
    for hit in self.hits:
      if hit.subject['revcomp']:
        raise FormatError('Invalid LAV: Subject sequence '+hit.subject['name']
          +'cannot be reverse complemented.')
      for alignment in hit.alignments:
        alignment = self._convert_segment(alignment, hit)
        for block in alignment.blocks:
          block = self._convert_segment(block, hit)
    self.converted = True


  def _convert_segment(self, segment, hit):
    """Convert any segment formatted like an LavAlignment or LavBlock.
    Specifically, anything with "subject" and "query" attributes which are dicts
    containing "begin" and "end" values. The "query" dict also must have a
    "revcomp" value."""
    subj_start = hit.subject['begin']
    quer_start = hit.query['begin']
    quer_end = hit.query['end']
    segment.subject['begin'] = segment.subject['begin'] + subj_start - 1
    segment.subject['end'] = segment.subject['end'] + subj_start - 1
    if hit.query['revcomp']:
      begin_temp = segment.query['begin']
      segment.query['begin'] = quer_end - segment.query['end'] + 1
      segment.query['end'] = quer_end - begin_temp + 1
    else:
      segment.query['begin'] = segment.query['begin'] + quer_start - 1
      segment.query['end'] = segment.query['end'] + quer_start - 1
    return segment


  def _stanza_start(self, line):
    """Detect if this line begins a stanza. If so, return the type (a single
    letter). If not, return None. The line should already be stripped."""
    fields = line.split()
    if len(fields) == 2 and fields[1] == '{':
      # stanza label is a single alphabetic character?
      if len(fields[0]) == 1 and 97 <= ord(fields[0].lower()) <= 122:
        return fields[0]


  def _parse_s(self, line, current_hit, stanza_line):
    fields = line.split()
    if len(fields) < 5:
      raise FormatError('Invalid LAV: LAV file must include all 5 fields '
        +'in "s" stanzas.')
    try:
      filename = fields[0].strip('"')
      begin    = int(fields[1])
      end      = int(fields[2])
      revcomp  = fields[3] == '1'
      seqnum   = int(fields[4])
    except ValueError:
      raise FormatError('Invalid LAV: Problem in "s" stanza: either a file'
        +'name with a space, or a non-integer coordinate.')
    if stanza_line == 1:
      current_hit.subject['filename'] = filename
      current_hit.subject['begin']    = begin
      current_hit.subject['end']      = end
      current_hit.subject['revcomp']  = revcomp
      current_hit.subject['seqnum']   = seqnum
    elif stanza_line == 2:
      current_hit.query['filename'] = filename
      current_hit.query['begin']    = begin
      current_hit.query['end']      = end
      current_hit.query['revcomp']  = revcomp
      current_hit.query['seqnum']   = seqnum
    return current_hit


  def _parse_h(self, line, current_hit, stanza_line):
    # remove quotes, '>' FASTA character, and reverse complement string
    name = line.strip('"')
    if name[0] == '>':
      name = name[1:]
    if name[len(name)-21:] == ' (reverse complement)':
      name = name[:len(name)-21]
    # get FASTA id
    fields = name.split()
    if fields:
      identifier = fields[0]
    else:
      identifier = ''
    if stanza_line == 1:
      current_hit.subject['id'] = identifier
      current_hit.subject['name'] = name
    elif stanza_line == 2:
      current_hit.query['id'] = identifier
      current_hit.query['name'] = name
    return current_hit


  def _parse_a(self, line, current_alignment, stanza_line):
    fields = line.split()
    if stanza_line == 1:
      if not (len(fields) == 2 and fields[0] == 's'):
        raise FormatError('Invalid LAV: Error in "s" line of "a" stanza.')
      current_alignment.score = int(fields[1])
    elif stanza_line == 2:
      if not (len(fields) == 3 and fields[0] == 'b'):
        raise FormatError('Invalid LAV: Error in "b" line of "a" stanza.')
      current_alignment.subject['begin'] = int(fields[1])
      current_alignment.query['begin']   = int(fields[2])
    elif stanza_line == 3:
      if not (len(fields) == 3 and fields[0] == 'e'):
        raise FormatError('Invalid LAV: Error in "e" line of "a" stanza.')
      current_alignment.subject['end'] = int(fields[1])
      current_alignment.query['end']   = int(fields[2])
    elif stanza_line >= 4:
      if not (len(fields) == 6 and fields[0] == 'l'):
        raise FormatError('Invalid LAV: Error in "l" line of "a" stanza.')
      block = LavBlock(current_alignment)
      block.subject['begin'] = int(fields[1])
      block.query['begin']   = int(fields[2])
      block.subject['end']   = int(fields[3])
      block.query['end']     = int(fields[4])
      block.identity         = int(fields[5])
      current_alignment.blocks.append(block)
    return current_alignment


class LavHit(object):

  def __init__(self, reader):
    self.parent = reader
    self.subject = {}
    self.query = {}
    self.alignments = []

  def __iter__(self):
    return iter(self.alignments)

  def __len__(self):
    return len(self.alignments)


class LavAlignment(object):

  def __init__(self, hit):
    self.parent = hit
    self.score = None
    self.subject = {}
    self.query = {}
    self.blocks = []

  def __iter__(self):
    return iter(self.blocks)

  def __len__(self):
    return len(self.blocks)


class LavBlock(object):

  def __init__(self, alignment):
    self.parent = alignment
    self.identity = None
    self.subject = {}
    self.query = {}

  def __len__(self):
    return self.subject['end'] - self.subject['begin'] + 1

