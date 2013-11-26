#!/usr/bin/env python

# TODO: Add a generator that reads an arbitrary number of characters at a time.
# Detect end of record with
#   header = seq.index('>')
#   if seq[header-1] == '\n':
#     break
#   elif header == 0 and lastchar = '\n':
#     break
# And then remove newlines with seq.replace('\n', '')
#
# TODO: Turn class into just a function
# See if this works:
#   def fasta_bases(filepath, which_seq=None, restart=True):
#     if restart:
#       filehhandle = open(filepath, 'rU')
#       reading = False
#       restart = False
# Then you'd call it like:
#   for base in fasta_bases('chrM.fa'):
#     print base
class FastaBaseGenerator:

  def __init__(self, filepath):
    self.filehandle = open(filepath, 'rU')
    self.reading = False

  def new(self, which_seq=None):

    seqnum = None
    seqid = None
    if isinstance(which_seq, int):
      seqnum = which_seq
    if isinstance(which_seq, str):
      seqid = which_seq

    # read until correct FASTA header
    record = 0
    while not self.reading:
      line_orig = self.filehandle.readline()
      if not line_orig:
        break
      line = line_orig.rstrip('\r\n')
      if not line:
        continue # allow blank lines
      if line[0] == '>':
        if seqnum:
          record+=1
          if record >= seqnum:
            self.reading = True
        elif seqid:
          record_id = line.split()[0]
          if record_id == '>'+seqid:
            self.reading = True
        else:
          self.reading = True

    newline = True
    while self.reading:
      base = self.filehandle.read(1)
      if not base:
        self.reading = False
      elif base == '\n':
        newline = True
      elif base == '>':
        if newline:
          self.reading = False
        else:
          newline = False
          yield base
      else:
        newline = False
        yield base


def main():
  import sys
  fbg = FastaBaseGenerator(sys.argv[1])
  base_generator = fbg.new()
  count = 0
  for base in base_generator:
    count+=1
    if count > 20:
      break
    print str(count)+": "+base

if __name__ == '__main__':
  main()