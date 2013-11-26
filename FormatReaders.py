#!/usr/bin/env python

class FormatException(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


# TODO:
#   Turn into iterator
#   Parse header in __init__()
#   Key sample fields with sample names
#   Delay parsing of extra fields until requested
class VCFReader(object):
  """A simple VCF parser which can read Naive Variant Caller output."""

  def __init__(self, filehandle):
    """Pass in a filehandle open in 'rU' mode (or at least 'r')"""
    self._filehandle = filehandle
    self._in_header = True
    self._header = ""
    self._sample_names = []

  def new(self):

    line_num = 0
    self._sample_names = []
    for line in self._filehandle:
      line_num+=1
      line = line.rstrip('\r\n')

      if self._in_header:
        if line[0] == '#':
          self._header += line+'\n'
          if line[0:6].upper() == '#CHROM':
            self._sample_names = line.split('\t')[9:]
          continue
        else:
          if self._sample_names:
            self._in_header = False
          else:
            raise FormatException("Invalid VCF: failed on line "+str(line_num))

      yield VCFPosition(line, line_num=line_num)

  def get_header(self):
    return self._header

  def get_sample_names(self):
    return self._sample_names

  def set_sample_names(self, sample_names):
    if isinstance(sample_names, list):
      self._sample_names = sample_names
      return True
    else:
      return False


class VCFPosition(object):

  def __init__(self, line, line_num=None):
    if line_num:
      self._line_num = str(line_num)
    else:
      self._line_num = "[NO NUMBER]"

    columns = line.split('\t')
    if len(columns) < 10:
      raise FormatException("Invalid VCF: too few columns in line "
        +self._line_num)

    self._chrom  = columns[0]
    try:
      self._pos  = int(columns[1])
    except ValueError:
      raise FormatException("Invalid VCF: non-integer POS in line "
        +self._line_num)
    self._id = columns[2]
    self._ref = columns[3]
    if columns[4] == '.':
      self._alt = []
    else:
      self._alt = columns[4].split(',')
    self._qual = columns[5]
    self._filter = columns[6]
    self._info = self._parse_info(columns[7])
    self._genotypes = self._parse_genotypes(columns[8], columns[9:])


  def _parse_info(self, info_string):
    info = {}

    for keyvalue in info_string.split(';'):
      try:
        (key, value) = keyvalue.split('=')
      except ValueError:
        raise FormatException("Invalid VCF: bad INFO field in line "
          +self._line_num)
      info[key] = value.split(',')

    return info


  def _parse_genotypes(self, format, samples):
    genotypes = []
    format_strings = format.split(':')
    
    for sample in samples:
      sample_strings = sample.split(':')
      if len(format_strings) != len(sample_strings):
        raise FormatException("Invalid VCF: FORMAT does not match SAMPLE "
          +"in line "+self._line_num)
      genotype = dict(zip(format_strings, sample_strings))
      genotypes.append(genotype)

    return genotypes


  def get_varcounts(self):
    varcounts = []

    for genotype in self._genotypes:
      varcount = {}

      try:
        varcount_strings = genotype['NC'].split(',')
      except KeyError:
        raise FormatException("Invalid VCF: may not be Naive Variant Counter "
          "output (line "+self._line_num+")")
      
      for varcount_string in varcount_strings:
        # the last one will always be empty
        if not varcount_string:
          continue
        vcfields = varcount_string.split('=')
        try:
          variant = vcfields[0]
          count = int(vcfields[1])
        except (IndexError, ValueError):
          raise FormatException("Invalid VCF: may not be Naive Variant Counter "
          "output (line "+self._line_num+")")
        varcount[variant] = count

      varcounts.append(varcount)

    return varcounts


  def get_line_num(self):
    return self._line_num

  def get_chrom(self):
    return self._chrom

  def get_pos(self):
    return self._pos

  def get_id(self):
    return self._id

  def get_ref(self):
    return self._ref

  def get_alt(self):
    return self._alt

  def get_qual(self):
    return self._qual

  def get_filter(self):
    return self._filter

  def get_info(self):
    return self._info

  def get_format(self):
    return self._format

  def get_genotypes(self):
    return self._genotypes

  def set_alt(self, alt):
    if isinstance(alt, list):
      self._alt = alt
      return True
    else:
      return False

  def set_genotypes(self, genotypes):
    if isinstance(genotypes, list):
      for genotype in genotypes:
        if not isinstance(genotype, dict):
          return False
      self._genotypes = genotypes
      return True
    else:
      return False


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
class FastaBaseGenerator(object):

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