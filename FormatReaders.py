#!/usr/bin/env python

class FormatException(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


#TODO: cache attributes instead of re-generating in getter
class VCFReader(object):
  """A simple VCF parser which can read Naive Variant Caller output."""

  def __iter__(self):
    return self

  def __init__(self, filehandle):
    """Pass in a filehandle open in 'rU' mode (or at least 'r')"""
    self._filehandle = filehandle
    self._in_header = True
    self._header = ""
    self._sample_names = []

    self._line_num = 0
    while self._in_header:
      line = self._filehandle.next()
      self._line_num+=1
      line = line.rstrip('\r\n')

      if line[0] == '#':
        self._header += line+'\n'
        if line[0:6].upper() == '#CHROM':
          self._sample_names = line.split('\t')[9:]
          self._in_header = False
      else:
        raise FormatException("Invalid VCF: invalid header (line "
          +self._line_num+")")


  def next(self):

    line = ""
    # allow empty lines
    while not line:
      line = self._filehandle.next()
      self._line_num+=1
      line = line.rstrip('\r\n')

    if self._in_header or line[0] == '#':
      raise FormatException("Invalid VCF: late header at line "+self._line_num)
    
    return VCFSite(line, self)


  def get_line_num(self):
    return self._line_num

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


class VCFSite(object):

  def __init__(self, line, reader):
    if reader.get_line_num():
      self._line_num = str(reader.get_line_num())
    else:
      self._line_num = "[NO NUMBER]"

    self._columns = line.split('\t')
    if len(self._columns) < 10:
      raise FormatException("Invalid VCF: too few columns in line "
        +self._line_num)

    self._reader = reader

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
    genotypes = {}
    format_strings = format.split(':')
    
    for (sample, sample_name) in zip(samples, self._reader.get_sample_names()):
      sample_strings = sample.split(':')
      if len(format_strings) != len(sample_strings):
        raise FormatException("Invalid VCF: FORMAT does not match SAMPLE "
          +"in line "+self._line_num)
      genotype = dict(zip(format_strings, sample_strings))
      genotypes[sample_name] = genotype

    return genotypes


  def _parse_varcounts(self, genotypes, stranded=True):
    varcounts = {}

    for sample_name in genotypes:
      varcount = {}

      try:
        varcount_strings = genotypes[sample_name]['NC'].split(',')
      except KeyError:
        raise FormatException("Invalid VCF: may not be Naive Variant Caller "
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
          raise FormatException("Invalid VCF: may not be Naive Variant Caller "
          "output (line "+self._line_num+")")
        if not stranded:
          variant = variant.lstrip('+-')
        varcount[variant] = count + varcount.get(variant, 0)

      varcounts[sample_name] = varcount

    return varcounts


  def _sum_coverages(self, varcounts):
    coverages = {}

    for sample_name in varcounts:
      total = 0
      varcount = varcounts[sample_name]
      for variant in varcount:
        total += varcount[variant]
      coverages[sample_name] = total

    return coverages


  def get_line_num(self):
    return self._line_num

  def get_chrom(self):
    if self._columns[0] == '.':
      return None
    else:
      return self._columns[0]

  def get_pos(self):
    try:
      return int(self._columns[1])
    except ValueError:
      raise FormatException("Invalid VCF: non-integer POS in line "
        +self._line_num)

  def get_id(self):
    if self._columns[2] == '.':
      return None
    else:
      return self._columns[2]

  def get_ref(self):
    if self._columns[3] == '.':
      return None
    else:
      return self._columns[3]

  def get_alt(self):
    if self._columns[4] == '.':
      return []
    else:
      return self._columns[4].split(',')

  def get_qual(self):
    if self._columns[5] == '.':
      return None
    else:
      return self._columns[5]

  def get_filter(self):
    if self._columns[6] == '.':
      return None
    else:
      return self._columns[6]

  def get_info(self):
    return self._parse_info(self._columns[7])

  def get_genotypes(self):
    return self._parse_genotypes(self._columns[8], self._columns[9:])

  def get_varcounts(self, stranded=True):
    return self._parse_varcounts(self.get_genotypes(), stranded=stranded)

  def get_coverages(self):
    return self._sum_coverages(self.get_varcounts())


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

  def set_varcounts(self, varcounts):
    if isinstance(varcounts, list):
      for varcount in varcounts:
        if not isinstance(varcount, dict):
          return False
      self._varcounts = varcounts
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
#
# TODO: check out https://pypi.python.org/pypi/pyfasta/
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