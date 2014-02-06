#!/usr/bin/env python
__version__ = '0.5'


class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


class FastaLineGenerator(object):
  """A simple FASTA parser that only reads a line at a time into memory.
  Usage:
  fasta = FastaLineGenerator('/home/user/sequence.fasta')
  for line in fasta:
    print "There is a sequence with this FASTA identifier: "+fasta.id
    print "It has a line with this sequence: "+line
  """

  def __init__(self, filepath):
    self.filehandle = open(filepath, 'rU')
    self.name = None
    self.id = None

  def __iter__(self):
    return self.new()

  def new(self):
    while True:
      line_raw = self.filehandle.readline()
      if not line_raw:
        raise StopIteration
      line = line_raw.strip()
      if not line:
        continue # allow empty lines
      if line[0] == '>':
        self.name = line[1:]  # remove ">"
        if self.name:
          self.id = self.name.split()[0]
        else:
          self.id = ''
        continue
      else:
        yield line


#TODO: see 0notes.txt
class FastaBaseGenerator(object):
  """For when you absolutely have to read one base at a time. VERY SLOW.
  Usage:
  fasta = FastaBaseGenerator('/home/user/sequence.fasta')
  for base in fasta:
    print "There is a sequence with this FASTA identifier: "+fasta.id
    print "This is the next base from it: "+base
  """

  def __init__(self, filepath):
    self.filehandle = open(filepath, 'rU')
    self.header = False
    self.name = None
    self.id = None
    self._in_id = None

  def __iter__(self):
    return self.new()

  def new(self):

    newline = True
    while True:
      base = self.filehandle.read(1)
      if not base:
        raise StopIteration
      elif base == '\n':
        newline = True
        self.header = False
      elif newline and base == '>':
        newline = False
        self.header = True
        self._in_id = True
        self.name = ''
        self.id = ''
      elif self.header:
        if self._in_id:
          if base.isspace():
            self._in_id = False
          else:
            self.id += base
        self.name += base
      else:
        newline = False
        yield base

