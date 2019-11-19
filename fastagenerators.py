#!/usr/bin/env python3
import collections.abc
import os
import pathlib
__version__ = '0.10'


class UsageError(Exception):
  pass


class Reader:
  """Base class for all other parsers."""
  def __init__(self, input, **kwargs):
    self.input = input
    self.input_type = detect_input_type(input)
    if self.input_type not in ('path', 'str', 'file', 'generator'):
      raise ValueError(f'Input object {input!r} not a file, path, string, or generator.')
    for key, value in kwargs.items():
      setattr(self, key, value)
  def get_input_iterator(self):
    if self.input_type == 'str':
      return open(self.input)
    elif self.input_type == 'path':
      return self.input.open()
    else:
      return self.input


class FastaSeqBuffered(object):
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


class FastaLineBuffered(Reader):
  """A simple FASTA parser that only reads a line at a time into memory.
  Usage:
  fasta = FastaLineBuffered('/home/user/sequence.fasta')
  for sequence in fasta:
    print(f'There is a sequence with this FASTA identifier: {sequence.id}')
    print(f'Its full name is {sequence.name}.')
    for line in sequence:
      print(f'It has a line with this sequence: {line}')
  All strings (the line, id, and name) are stripped, and should not end in a
  newline.
  """

  def __init__(self, input):
    super().__init__(input)
    if self.input_type is 'str' and not os.path.isfile(self.input):
      raise IOError(f'File not found: "{self.input!r}"')
    self._sequence = None
    self._initialized = False

  def __iter__(self):
    self._input_iterator = self.get_input_iterator()
    # Prime it by reading until the first sequence header.
    line = None
    while line is None or not line.startswith('>'):
      try:
        line_raw = next(self._input_iterator)
      except StopIteration:
        self._cleanup()
        break
      line = line_raw.strip()
    self._last_line = line
    self._initialized = True
    return self

  def __next__(self):
    if not self._initialized:
      raise UsageError('Iterator not initialized!')
    if self._sequence is not None:
      if not self._sequence.done:
        # If the caller requested the next sequence before the last was done generating lines,
        # finish up the previous one by reading lines until we find the next sequence.
        self._sequence.fast_forward()
      next_line = self._sequence._last_line
    else:
      next_line = self._last_line
    if next_line is None:
      self._cleanup()
      raise StopIteration
    seq_id, seq_name = parse_header(next_line)
    try:
      self._sequence = SequenceLineBuffered(seq_id, seq_name, self._input_iterator)
    except StopIteration:
      self._cleanup()
      raise
    return self._sequence

  def _cleanup(self):
    if self.input_type in ('path', 'file', 'str'):
      self._input_iterator.close()

  def extract(self, start, end, chrom=None):
    """Extract a subsequence based on a start and end coordinate.
    If chrom is not supplied, it will default to the first sequence (record)
    encountered in the FASTA file.
    See SequenceLineBuffered.extract() for details on the rest."""
    for sequence in self.sequences():
      if chrom is not None and sequence.id != chrom:
        continue
      return sequence.extract(start, end)


class SequenceLineBuffered:

  def __init__(self, seq_id, seq_name, input_iterator):
    self.name = seq_name
    self.id = seq_id
    self.done = False
    self._input_iterator = input_iterator
    self._last_line = None

  def __iter__(self):
    return self

  def __next__(self):
    if self.done:
      raise UsageError(f'Tried to read from a finished {type(self).__name__}.')
    line = self._get_next_line()
    if line.startswith('>'):
      self._last_line = line
      self.done = True
      raise StopIteration
    else:
      return line

  def _get_next_line(self):
    """Return the next nonempty line or raise StopIteration if there are none."""
    line = None
    while not line:
      try:
        line_raw = next(self._input_iterator)
      except StopIteration:
        self.done = True
        raise StopIteration
      line = line_raw.strip()
    return line

  def fast_forward(self):
    while True:
      try:
        next(self)
      except StopIteration:
        break

  def bases(self):
    """Generator that yields single bases, while still reading a whole line at
    a time underneath.
    This should be the best of both worlds: it yields a base at a time, but it
    reads a line at a time from the file so it's not slow as molasses."""
    for line in self:
      for base in line:
        yield base

  def extract(self, start, end):
    """Extract a subsequence based on a start and end coordinate.
    The start and end are inclusive, 1-based. If the end coordinate is beyond the
    end of the chromosome, the returned sequence will be truncated to the end of the
    chromosome. If the start coordinate is beyond the end of the chromosome, an
    empty string will be returned."""
    outseq = ''
    line_start = 1
    for line in self:
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
class FastaBaseBuffered(object):
  """For when you absolutely have to read one base at a time (e.g. very long lines). VERY SLOW.
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


def detect_input_type(input):
  if isinstance(input, pathlib.Path):
    return 'path'
  elif isinstance(input, str):
    return 'str'
  elif isinstance(input, collections.abc.Generator):
    return 'generator'
  elif hasattr(input, 'read') and hasattr(input, 'close'):
    return 'file'
  else:
    return None


def parse_header(line):
  seq_id = None
  seq_name = None
  if line.startswith('>'):
    seq_name = line[1:]
    if seq_name:
      seq_id = seq_name.split()[0]
    else:
      seq_id = None
  return seq_id, seq_name
