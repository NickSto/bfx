import os
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


def getparser(input, filetype='fasta'):
  # Detect whether the input is an open file or a path.
  # Return the appropriate reader.
  if filetype == 'fasta':
    return FastaReader(input)
  elif filetype == 'fastq':
    return FastqReader(input)
  elif filetype == 'sam':
    return SamReader(input)
  elif filetype == 'tsv':
    return TsvReader(input)
  else:
    raise ValueError('Unrecognized format: {!r}'.format(filetype))


def detect_input_type(obj):
  """Is this an open filehandle, or is it a file path (string)?"""
  try:
    os.path.isfile(obj)
    return 'path'
  except TypeError:
    return 'file'


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
  def __init__(self, input):
    self.input = input
    self.input_type = detect_input_type(input)
    if self.input_type not in ('path', 'file'):
      raise ValueError('Input object {!r} not file-like or string-like.'.format(input))
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
    if self.input_type == 'path':
      filehandle = open(self.input)
    elif self.input_type == 'file':
      filehandle = self.input
    try:
      for line in filehandle:
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
    finally:
      if self.input_type == 'path':
        filehandle.close()


class SamReader(Reader):
  """A simple SAM parser.
  Assumptions:
  Lines starting with "@" with 3 fields are headers. All others are alignments.
  All alignment lines have 11 or more fields. Other lines will be skipped.
  """
  def parser(self):
    if self.input_type == 'path':
      filehandle = open(self.input)
    elif self.input_type == 'file':
      filehandle = self.input
    try:
      for line in filehandle:
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
    finally:
      if self.input_type == 'path':
        filehandle.close()


class FastaReader(Reader):
  """A simple FASTA parser that reads one sequence at a time into memory."""
  def parser(self):
    if self.input_type == 'path':
      filehandle = open(self.input)
    elif self.input_type == 'file':
      filehandle = self.input
    try:
      read = Read()
      while True:
        line_raw = filehandle.readline()
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
    finally:
      if self.input_type == 'path':
        filehandle.close()


class FastqReader(Reader):
  """A simple FASTQ parser. Can handle multi-line sequences, though."""
  def parser(self):
    if self.input_type == 'path':
      filehandle = open(self.input)
    elif self.input_type == 'file':
      filehandle = self.input
    try:
      read = Read()
      state = 'header'
      while True:
        line_raw = filehandle.readline()
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
    finally:
      if self.input_type == 'path':
        filehandle.close()
