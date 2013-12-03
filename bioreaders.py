#!/usr/bin/env python
# requires Python 2.7
from collections import OrderedDict

class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)


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
        raise FormatError("Invalid VCF: invalid header (line "
          +str(self._line_num)+")")


  def next(self):

    line = ""
    # allow empty lines
    while not line:
      line = self._filehandle.next()
      self._line_num+=1
      line = line.rstrip('\r\n')

    if self._in_header or line[0] == '#':
      raise FormatError("Invalid VCF: late header at line "
        +str(self._line_num))
    
    site = VCFSite()
    site.parse_line(line)
    site.set_line_num(self._line_num)
    site.set_sample_names(self._sample_names)
    return site


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

  def __init__(self):
    # internal use
    self._initialized = False
    self._sample_names = None
    self._line_num = '[N/A]'
    self._modified = [False] * 10
    # data
    self._columns = None
    self._chrom = None
    self._pos = None
    self._id = None
    self._ref = None
    self._alt = None
    self._qual = None
    self._filter = None
    self._info = None
    self._genotypes = None
    self._varcounts_stranded = None
    self._varcounts_unstranded = None
    self._variants_stranded = None
    self._variants_unstranded = None
    self._coverages = None


  def parse_line(self, line):
    self._columns = line.split('\t')
    if len(self._columns) < 10:
      raise FormatError("Invalid VCF: too few columns in line "+self._line_num)
    self._initialized = True


  def get_chrom(self):
    if self._chrom is None and self._initialized:
      if self._columns[0] == '.':
        self._chrom = None
      else:
        self._chrom = self._columns[0]
    return self._chrom

  def get_pos(self):
    if self._pos is None and self._initialized:
      try:
        self._pos = int(self._columns[1])
      except ValueError:
        raise FormatError("Invalid VCF: non-integer POS in line "
          +self._line_num)
    return self._pos

  def get_id(self):
    if self._id is None and self._initialized:
      if self._columns[2] == '.':
        self._id = None
      else:
        self._id = self._columns[2]
    return self._id

  def get_ref(self):
    if self._ref is None and self._initialized:
      if self._columns[3] == '.':
        self._ref = None
      else:
        self._ref = self._columns[3]
    return self._ref

  def get_alt(self):
    if self._alt is None and self._initialized:
      if self._columns[4] == '.':
        self._alt = []
      else:
        self._alt = self._columns[4].split(',')
    return self._alt

  def get_qual(self):
    if self._qual is None and self._initialized:
      if self._columns[5] == '.':
        self._qual = None
      else:
        try:
          self._qual = int(self._columns[5])
        except ValueError:
          self._qual = float(self._columns[5])
        except ValueError:
          raise FormatError("Invalid VCF: non-numeric QUAL in line "
            +self._line_num)
    return self._qual

  def get_filter(self):
    if self._filter is None and self._initialized:
      if self._columns[6] == '.':
        self._filter = None
      elif self._columns[6].upper() == 'PASS':
        self._filter = True
      else:
        self._filter = self._columns[6].split(';')
    return self._filter

  def get_info(self):
    if self._info is None and self._initialized:
      self._info = self._parse_info(self._columns[7])
    return self._info

  def get_genotypes(self):
    if self._genotypes is None and self._initialized:
      self._genotypes = self._parse_genotypes(self._columns[8],
        self._columns[9:])
    return self._genotypes

  def get_varcounts(self, stranded=False):
    if stranded:
      if self._varcounts_stranded is None and self._initialized:
        self._varcounts_stranded = self._parse_varcounts(self.get_genotypes(),
          stranded=True)
      return self._varcounts_stranded
    else:
      if self._varcounts_unstranded is None and self._initialized:
        self._varcounts_unstranded = self._parse_varcounts(self.get_genotypes(),
          stranded=False)
      return self._varcounts_unstranded

  def get_variants(self, stranded=False):
    if stranded:
      if self._variants_stranded is None and self._initialized:
        self._variants_stranded = self._variants_list(self.get_varcounts(
          stranded=True))
      return self._variants_stranded
    else:
      if self._variants_unstranded is None and self._initialized:
        self._variants_unstranded = self._variants_list(self.get_varcounts(
          stranded=False))
      return self._variants_unstranded

  def get_coverages(self):
    if self._coverages is None and self._initialized:
      self._coverages = self._sum_coverages(self.get_varcounts())
    return self._coverages

  def get_sample_names(self):
    return self._sample_names

  def get_line_num(self):
    return self._line_num


  def set_chrom(self, chrom):
    if isinstance(chrom, basestring) or chrom is None:
      self._chrom = chrom
      self._modified[0] = True
      return True
    else:
      return False

  def set_pos(self, pos):
    if isinstance(pos, int):
      self._pos = pos
      self._modified[1] = True
      return True
    else:
      return False

  def set_id(self, id):
    if isinstance(id, basestring) or id is None:
      self._id = id
      self._modified[2] = True
      return True
    else:
      return False

  def set_ref(self, ref):
    if isinstance(ref, basestring) or ref is None:
      self._ref = ref
      self._modified[3] = True
      return True
    else:
      return False

  #TODO: Alter the NVC data in the INFO and genotypes to be consistent with new
  #  list of ALTs
  def set_alt(self, alt):
    if isinstance(alt, list) or alt is None:
      self._alt = alt
      self._modified[4] = True
      return True
    else:
      return False

  def set_qual(self, qual):
    if isinstance(qual, (int, long, float)) or qual is None:
      self._qual = qual
      self._modified[5] = True
      return True
    else:
      return False

  def set_filter(self, filter):
    if isinstance(filter, list) or filter is None or filter is True:
      self._filter = filter
      self._modified[6] = True
      return True
    else:
      return False

  def set_info(self, info):
    if isinstance(info, dict) or info is None:
      self._info = info
      self._modified[7] = True
      return True
    else:
      return False

  def set_genotypes(self, genotypes):
    if isinstance(genotypes, dict):
      for genotype in genotypes.values():
        if not isinstance(genotype, dict):
          return False
      self._genotypes = genotypes
      self._modified[8] = True
      return True
    else:
      return False

  #TODO: implement (remember to roll varcounts changes into genotypes)
  # def set_varcounts(self, varcounts):

  def set_sample_names(self, sample_names):
    if isinstance(sample_names, list):
      self._sample_names = sample_names
      return True
    else:
      return False

  def set_line_num(self, line_num):
    self._line_num = str(line_num)
    return True


  def variant_to_alt(self, variant):
    variant = variant.lstrip('+-')
    ref = self.get_ref()
    ref_tail = ref[1:]
    if variant[0] == 'd':
      try:
        delength = int(variant[1:])
      except ValueError:
        raise FormatError("Invalid variant string")
      return ref[0]+ref_tail[delength:]
    else:
      return variant+ref_tail


  def split():
    """Split a multi-sample VCFSite into multiple VCFSites, one per sample."""

    sites = []
    for sample_name in self._sample_names:
      sites.append(sample_name)


  def __str__(self):
    """Returns a VCF line with the site's current data (no newline)"""
    if True not in self._modified:
      return '\t'.join(self._columns)

    line = []
    line.append(self.get_chrom())
    line.append(self.get_pos())
    line.append(self.get_id())
    line.append(self.get_ref())
    line.append(','.join(self.get_alt()))
    line.append(self.get_qual())
    # FILTER column
    if self._modified[6]:
      filter = self.get_filter()
      if filter is True:
        line.append("PASS")
      else:
        line.append(';'.join([str(val) for val in filter]))
    else:
      line.append(self._columns[6])
    # INFO column
    if self._modified[7]:
      info = self.get_info()
      fields = []
      for key in info:
        data_string = ','.join([str(val) for val in info[key]])
        fields.append(str(key)+'='+data_string)
      line.append(';'.join(fields))
    else:
      line.append(self._columns[7])
    # FORMAT and sample columns
    if self._modified[8]:
      genotypes = self.get_genotypes()
      data = genotypes.values()
      line.append(':'.join(data[0].keys()))  # FORMAT
      for sample_name in genotypes:
        sample_column = ':'.join(genotypes[sample_name].values())
        line.append(sample_column)
    else:
      line.extend(self._columns[8:])

    for i in range(len(line)):
      if line[i] is None:
        line[i] = '.'
    return '\t'.join([str(field) for field in line])


  def _parse_info(self, info_string):
    info = OrderedDict()

    for keyvalue in info_string.split(';'):
      if '=' in keyvalue:
        (key, value) = keyvalue.split('=')
        info[key] = value.split(',')
      else:
        info[keyvalue] = True

    return info


  def _parse_genotypes(self, format, samples):
    genotypes = OrderedDict()
    format_strings = format.split(':')
    
    for (sample, sample_name) in zip(samples, self._sample_names):
      sample_strings = sample.split(':')
      if len(format_strings) != len(sample_strings):
        raise FormatError("Invalid VCF: FORMAT does not match SAMPLE in line "
          +self._line_num)
      genotype = OrderedDict(zip(format_strings, sample_strings))
      genotypes[sample_name] = genotype

    return genotypes


  def _parse_varcounts(self, genotypes, stranded=False):
    varcounts = OrderedDict()

    for sample_name in genotypes:
      varcount = OrderedDict()

      try:
        varcount_strings = genotypes[sample_name]['NC'].split(',')
      except KeyError:
        raise FormatError("Invalid VCF: may not be Naive Variant Caller "
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
          raise FormatError("Invalid VCF: may not be Naive Variant Caller "
          "output (line "+self._line_num+")")
        if not stranded:
          variant = variant.lstrip('+-')
        varcount[variant] = count + varcount.get(variant, 0)

      varcounts[sample_name] = varcount

    return varcounts


  def _variants_list(self, varcounts):
    variants = OrderedDict()
    for sample_name in varcounts:
      for variant in varcounts[sample_name]:
        variants[variant] = True
    return variants.keys()


  def _sum_coverages(self, varcounts):
    coverages = OrderedDict()

    for sample_name in varcounts:
      total = 0
      varcount = varcounts[sample_name]
      for variant in varcount:
        total += varcount[variant]
      coverages[sample_name] = total

    return coverages


#TODO: see 0notes.txt
class FastaBaseGenerator(object):

  def __init__(self, filepath):
    self.filehandle = open(filepath, 'rU')
    self.reading = False

  def new(self, which_seq=None):

    seqnum = None
    seqid = None
    if isinstance(which_seq, int):
      seqnum = which_seq
    if isinstance(which_seq, basestring):
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