"""A simple parser for FASTA, FASTQ, SAM, etc. Create generators that just return the read name and
sequence.
All format parsers follow this API:
  with open('sequence.fasta') as fasta:
    for read in getreads.getparser(fasta, filetype='fasta'):
      print "There is a sequence with this FASTA identifier: "+read.id
      print "Its sequence is "+read.seq
The properties of Read are:
  name: The entire FASTA header line, SAM column 1, etc.
  id:   The first whitespace-delimited part of the name.
  seq:  The sequence.
  qual: The quality scores (unless the format is FASTA).
"""


def getparser(filehandle, filetype='fasta'):
  if filetype == 'fasta':
    return FastaReader(filehandle)
  elif filetype == 'fastq':
    return FastqReader(filehandle)
  elif filetype == 'sam':
    return SamReader(filehandle)
  elif filetype == 'tsv':
    return TsvReader(filehandle)
  else:
    raise ValueError('Illegal argument: filetype="{}"'.format(filetype))


class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


class Read(object):
  def __init__(self, name='', seq='', id_='', qual=''):
    self.name = name
    self.seq = seq
    self.id = id_
    self.qual = qual


class Reader(object):
  """Base class for all other parsers."""
  def __init__(self, filehandle):
    self.filehandle = filehandle
  def __iter__(self):
    return self.parser()
  def bases(self):
    for read in self.parser():
      for base in read.seq:
        yield base


class TsvReader(Reader):
  """A parser for a simple tab-delimited format.
  Column 1: name
  Column 2: sequence
  Column 3: quality scores (optional)"""
  def parser(self):
    for line in self.filehandle:
      fields = line.rstrip('\r\n').split('\t')
      if len(fields) < 2:
        continue
      read = Read()
      read.name = fields[0]
      if read.name:
        read.id = read.name.split()[0]
      read.seq = fields[1]
      if len(fields) >= 3:
        read.qual = fields[2]
      yield read


class SamReader(Reader):
  """A simple SAM parser.
  Assumptions:
  Lines starting with "@" with 3 fields are headers. All others are alignments.
  All alignment lines have 11 or more fields. Other lines will be skipped.
  """
  def parser(self):
    for line in self.filehandle:
      fields = line.split('\t')
      if len(fields) < 11:
        continue
      # Skip headers.
      if fields[0].startswith('@') and len(fields[0]) == 3:
        continue
      read = Read()
      read.name = fields[0]
      if read.name:
        read.id = read.name.split()[0]
      read.seq = fields[9]
      read.qual = fields[10].rstrip('\r\n')
      yield read


class FastaReader(Reader):
  """A simple FASTA parser that reads one sequence at a time into memory."""
  def parser(self):
    read = Read()
    while True:
      line_raw = self.filehandle.readline()
      if not line_raw:
        if read.seq:
          yield read
        return
      line = line_raw.strip()
      # Allow empty lines.
      if not line:
        continue
      if line.startswith('>'):
        if read.seq:
          yield read
        read = Read()
        read.name = line[1:]  # remove ">"
        if read.name:
          read.id = read.name.split()[0]
        continue
      else:
        read.seq += line


class FastqReader(Reader):
  """A simple FASTQ parser. Can handle multi-line sequences, though."""
  def parser(self):
    read = Read()
    state = 'header'
    while True:
      line_raw = self.filehandle.readline()
      if not line_raw:
        if read.seq:
          yield read
        return
      line = line_raw.strip()
      # Allow empty lines.
      if not line:
        continue
      if state == 'header':
        if not line.startswith('@'):
          raise FormatError('line state = "header" but line does not start with "@"')
        if read.seq:
          yield read
        read = Read()
        read.name = line[1:]  # remove '@'
        if read.name:
          read.id = read.name.split()[0]
        state = 'sequence'
      elif state == 'sequence':
        if line.startswith('+'):
          state = 'plus'
        else:
          read.seq += line
      elif state == 'plus' or state == 'quality':
        state = 'quality'
        togo = len(read.seq) - len(read.qual)
        read.qual += line[:togo]
        # The end of the quality lines is when we have a quality string as long as the sequence.
        if len(read.qual) >= len(read.seq):
          state = 'header'
