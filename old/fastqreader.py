#!/usr/bin/env python
import os
__version__ = '0.5'


class FastqReadGenerator(object):
  """A simple FASTQ parser that returns reads one at a time.
  Handles multi-line read/quality values.
  Usage:
  fastq = FastqReadGenerator('/home/user/sequence.fq')
  for read in fastq:
    print "There is a read with this identifier: "+read.id
    print "(Its full name is "+read.name+".)"
    print "Its sequence is: "+read.seq
    print "Its quality is:  "+read.qual
  All values (id, name, seq, qual) are whitespace-stripped.
  """

  def __init__(self, filepath):
    if not os.path.isfile(filepath):
      raise IOError('File not found: "'+filepath+'"')
    self.filepath = filepath
    self.name = None
    self.id = None

  def __iter__(self):
    return self.reads()

  def reads(self):
    with open(self.filepath, 'rU') as filehandle:
      read = None
      line_type = 'first'
      for line_raw in filehandle:
        line = line_raw.strip()
        if not line:
          continue  # allow empty lines
        # Determine what kind of line we're in
        if line.startswith('@'):
          if line_type == 'first':
            line_type = 'name'
          elif line_type == 'plus':
            line_type = 'qual'
          elif line_type == 'qual':
            # Determine if it's another qual line or a name line.
            # If the quality scores observed so far already cover the whole read, we've seen all
            # the quality information already. It should be a name line.
            if len(read.qual) >= len(read.seq):
              line_type = 'name'
            else:
              line_type = 'qual'
          else:
            raise FormatError('"@" starts line in wrong context:\n'+line_raw)
        elif line.startswith('+'):
          if line_type == 'seq':
            line_type = 'plus'
          elif line_type == 'qual':
            pass
          else:
            raise FormatError('"+" starts line in wrong context:\n'+line_raw)
        elif line_type == 'name':
          line_type = 'seq'
        elif line_type == 'plus':
          line_type = 'qual'
        elif line_type == 'first':
          raise FormatError('First line must start with a "@":\n'+line_raw)
        else:
          raise FormatError('Invalid parser state: line_type "{}", first char "{}":\n{}'
                            .format(line_type, line[0], line_raw))
        if line_type == 'name':
          # Return the previous read.
          if read is not None:
            yield read
          read = Read()
          read.name = line[1:]  # remove ">"
          if read.name:
            read.id = read.name.split()[0]
          else:
            read.id = ''
        elif line_type == 'seq':
          read.seq += line
        elif line_type == 'qual':
          read.qual += line
      # Return the last read.
      if read is not None:
        yield read


class Read(object):
  def __init__(self):
    self.seq = ''
    self.qual = ''
    self.id = ''
    self.name = ''


class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)
