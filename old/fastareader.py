#!/usr/bin/env python
import os
__version__ = '0.9'


class FastaReadGenerator(object):
  """Read FASTA files and return one whole sequence at a time."""

  def __init__(self, filepath):
    self.line_generator = FastaLineGenerator(filepath)

  def __iter__(self):
    return self.reads()

  def reads(self):
    read = Read()
    read.name = None
    for line in self.line_generator:
      if self.line_generator.name != read.name:
        if read.name is not None:
          yield read
        read = Read()
        read.name = self.line_generator.name
        read.id = self.line_generator.id
      read.seq += line
    if read.name is not None:
      yield read


class FastaLineGenerator(object):
  """A simple FASTA parser that only reads a line at a time into memory.
  Usage:
  fasta = FastaLineGenerator('/home/user/sequence.fasta')
  for line in fasta:
    print "There is a sequence with this FASTA identifier: "+fasta.id
    print "(Its full name is "+fasta.name+".)"
    print "It has a line with this sequence: "+line
  All strings (the line, id, and name) are stripped, and should not end in a
  newline.
  """

  def __init__(self, filepath):
    if not os.path.isfile(filepath):
      raise IOError('File not found: "'+filepath+'"')
    self.filepath = filepath
    self.name = None
    self.id = None

  def __iter__(self):
    return self.lines()

  #TODO: Give some signal that we just finished a sequence. Otherwise, we can't validate that there
  #      aren't sequences with identical names one after another.
  def lines(self):
    with open(self.filepath, 'rU') as filehandle:
      for line_raw in filehandle:
        line = line_raw.strip()
        if not line:
          continue  # allow empty lines
        if line.startswith('>'):
          self.name = line[1:]  # remove ">"
          if self.name:
            self.id = self.name.split()[0]
          else:
            self.id = ''
          continue
        else:
          yield line

  def bases(self):
    """Generator that yields single bases, while still reading a whole line at
    a time underneath.
    This should be the best of both worlds: it yields a base at a time, but it
    reads a line at a time from the file so it's not slow as molasses."""
    for line in self.lines():
      for base in line:
        yield base

  def extract(self, start, end, chrom=None):
    """Extract a subsequence based on a start and end coordinate.
    The start and end are inclusive, 1-based. If chrom is not supplied, it will
    default to the first chromosome (record) encountered in the FASTA file.
    If the end coordinate is beyond the end of the chromosome, the returned
    sequence will be truncated to the end of the chromosome. If the start
    coordinate is beyond the end of the chromosome, an empty string will be
    returned."""
    outseq = ''
    line_start = 1
    for line in self:
      if chrom is not None and self.id != chrom:
        continue
      line_end = line_start + len(line) - 1
      # if we haven't encountered the start yet, keep searching
      if line_end < start:
        line_start = line_end + 1
        continue
      slice_start = max(start, line_start) - line_start
      slice_end = min(end, line_end) - line_start + 1
      outseq += line[slice_start:slice_end]
      # done? (on the last line?)
      if line_end >= end:
        break
      line_start = line_end + 1
    return outseq


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

class Read(object):
  def __init__(self):
    self.seq = ''
    self.id = ''
    self.name = ''
