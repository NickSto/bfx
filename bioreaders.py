#!/usr/bin/env python
# requires Python 2.7
__version__ = 'c0c4122'
from collections import OrderedDict
import copy

class FormatError(Exception):
  def __init__(self, message=None):
    if message:
      Exception.__init__(self, message)




class LavReader(object):
  """Parse an LAV file and provide an API for querying its fields.
  By default, the data will be represented exactly as it appears in the file
  (but with numbers as int types). To convert coordinates to a more intuitive
  system, call the convert() method.
    Data structures:
  LavReader.hits        A list of LavHits
  LavHit.subject        A dict containing data on the hit's subject sequence
  LavHit.query          A dict containing data on the hit's query sequence
    keys for subject and query:
    filename, id, name, seqnum, revcomp, begin, end 
  LavHit.alignments     A list of LavAlignments
  LavAlignment.query    A dict containing the query start and end coordinates
  LavAlignment.subject  A dict containing the subject start and end coordinates
    keys for subject and query:
    begin, end
  LavAlignment.score    An int: the score of the hit
  LavAlignment.blocks   A list of dicts, one for each gap-free block in the hit
    keys for blocks:
    subject_end, subject_begin, query_begin, query_end, length, identity
  """
  # Format assumptions:
  # This was written only to parse the output of LASTZ version 1.02.00.
  # Stanza starts and ends are on their own lines
  # - E.g. there will be nothing (except whitespace) before or after "a {" on
  #   the line in which it appears. The same goes for "}".
  # Stanza labels are single alphabetic characters
  # Sequence file names do not contain whitespace
  # h stanzas are present
  # s stanzas:
  # - rev_comp_flag's and sequence_number's are given
  # - only query sequences can be reverse complemented

  def __init__(self, filepath):
    self.hits = []
    
    stanza = ''
    stanza_done = True
    stanza_line = 0
    line_num = 0
    current_hit = LavHit()
    current_alignment = LavAlignment()
    with open(filepath, 'rU') as filehandle:
      for raw_line in filehandle:
        line = raw_line.strip()
        line_num+=1
        if not line:
          continue
        # at start of another stanza?
        if stanza_done:
          stanza = self._stanza_start(line)
          if stanza:
            stanza_done = False
            stanza_line = 0
          continue
        # at end of a stanza?
        if line == '}':
          # if end of an alignment block, add it to the current hit
          if stanza == 'a' and len(current_alignment) > 0:
            current_hit.alignments.append(current_alignment)
            current_alignment = LavAlignment()
          stanza_done = True
          stanza = ''
          continue
        # parse the stanzas
        stanza_line+=1
        if stanza == 's':
          # add previous hit to the list, start a new one
          if len(current_hit) > 0:
            self.hits.append(current_hit)
            current_hit = LavHit()
          # file names, input sequence info (revcomp, etc)
          try:
            current_hit = self._parse_s(line, current_hit, stanza_line)
          except FormatError as fe:
            fe.args = (fe.args[0]+' (line '+str(line_num)+')',)
            raise fe
        elif stanza == 'h':
          # sequence names
          current_hit = self._parse_h(line, current_hit, stanza_line)
        elif stanza == 'a':
          try:
            current_alignment = self._parse_a(line, current_alignment, stanza_line)
          except FormatError as fe:
            fe.args = (fe.args[0]+' (line '+str(line_num)+')',)
            raise fe
          except ValueError:
            raise FormatError('Invalid LAV: Non-integer encountered in "a" '
              +'stanza on line '+str(line_num)+'.')
    if len(current_alignment) > 0:
      current_hit.alignments.append(current_alignment)
    if len(current_hit) > 0:
      self.hits.append(current_hit)


  def convert(self):
    """Convert all coordinates to forward strand, 1-based, from the start.
    There cannot be any reverse-complemented subject sequences, or this will
    raise a FormatError."""
    for hit in self.hits:
      sstart = hit.subject['begin']
      qstart = hit.query['begin']
      qend = hit.query['end']
      if hit.subject['revcomp']:
        raise FormatError('Invalid LAV: Subject sequence '+hit.subject['name']
          +'cannot be reverse complemented.')
      for aln in hit.alignments:
        aln.subject['begin'] = aln.subject['begin'] + sstart - 1
        aln.subject['end'] = aln.subject['end'] + sstart - 1
        if hit.query['revcomp']:
          begin_temp = aln.query['begin']
          aln.query['begin'] = qend - aln.query['end'] + 1
          aln.query['end'] = qend - begin_temp + 1
        else:
          aln.query['begin'] = aln.query['begin'] + qstart - 1
          aln.query['end'] = aln.query['end'] + qstart - 1


  def _stanza_start(self, line):
    """Detect if this line begins a stanza. If so, return the type (a single
    letter). If not, return None. The line should already be stripped."""
    fields = line.split()
    if len(fields) == 2 and fields[1] == '{':
      # stanza label is a single alphabetic character?
      if len(fields[0]) == 1 and 97 <= ord(fields[0].lower()) <= 122:
        return fields[0]


  def _parse_s(self, line, current_hit, stanza_line):
    fields = line.split()
    if len(fields) < 5:
      raise FormatError('Invalid LAV: LAV file must include all 5 fields '
        +'in "s" stanzas.')
    try:
      filename = fields[0].strip('"')
      begin = int(fields[1])
      end = int(fields[2])
      revcomp = fields[3] == '1'
      seqnum = int(fields[4])
    except ValueError:
      raise FormatError('Invalid LAV: Problem in "s" stanza: either a file'
        +'name with a space, or a non-integer coordinate.')
    if stanza_line == 1:
      current_hit.subject['filename'] = filename
      current_hit.subject['begin']    = begin
      current_hit.subject['end']      = end
      current_hit.subject['revcomp']  = revcomp
      current_hit.subject['seqnum']   = seqnum
    elif stanza_line == 2:
      current_hit.query['filename'] = filename
      current_hit.query['begin']    = begin
      current_hit.query['end']      = end
      current_hit.query['revcomp']  = revcomp
      current_hit.query['seqnum']   = seqnum
    return current_hit


  def _parse_h(self, line, current_hit, stanza_line):
    # remove quotes, '>' FASTA character, and reverse complement string
    name = line.strip('"')
    if name[0] == '>':
      name = name[1:]
    if name[len(name)-21:] == ' (reverse complement)':
      name = name[:len(name)-21]
    # get FASTA id
    fields = name.split()
    if fields:
      identifier = fields[0]
    else:
      identifier = ''
    if stanza_line == 1:
      current_hit.subject['id'] = identifier
      current_hit.subject['name'] = name
    elif stanza_line == 2:
      current_hit.query['id'] = identifier
      current_hit.query['name'] = name
    return current_hit


  def _parse_a(self, line, current_alignment, stanza_line):
    fields = line.split()
    if stanza_line == 1:
      if not (len(fields) == 2 and fields[0] == 's'):
        raise FormatError('Invalid LAV: Error in "s" line of "a" stanza.')
      current_alignment.score = int(fields[1])
    elif stanza_line == 2:
      if not (len(fields) == 3 and fields[0] == 'b'):
        raise FormatError('Invalid LAV: Error in "b" line of "a" stanza.')
      current_alignment.subject['begin'] = int(fields[1])
      current_alignment.query['begin']   = int(fields[2])
    elif stanza_line == 3:
      if not (len(fields) == 3 and fields[0] == 'e'):
        raise FormatError('Invalid LAV: Error in "e" line of "a" stanza.')
      current_alignment.subject['end'] = int(fields[1])
      current_alignment.query['end']   = int(fields[2])
    elif stanza_line >= 4:
      if not (len(fields) == 6 and fields[0] == 'l'):
        raise FormatError('Invalid LAV: Error in "l" line of "a" stanza.')
      block = {}
      block['subject_begin'] = int(fields[1])
      block['query_begin']   = int(fields[2])
      block['subject_end']   = int(fields[3])
      block['query_end']     = int(fields[4])
      block['identity']      = int(fields[5])
      block['length']        = int(fields[3]) - int(fields[1]) + 1
      current_alignment.blocks.append(block)
    return current_alignment


class LavHit(object):

  def __init__(self):
    self.query = {}
    self.subject = {}
    self.alignments = []

  def __len__(self):
    return len(self.alignments)


class LavAlignment(object):

  def __init__(self):
    self.score = None
    self.query = {}
    self.subject = {}
    self.blocks = []

  def __len__(self):
    return len(self.blocks)



class VCFReader(object):
  """A simple VCF parser which can read Naive Variant Caller output."""

  def __iter__(self):
    return self

  def __init__(self, filehandle):
    """Pass in a filehandle open in 'rU' mode (or at least 'r')"""
    self._filehandle = filehandle
    self._in_header = True
    self._meta_header = ""
    self._column_header = ""
    self._sample_names = []

    self._line_num = 0
    while self._in_header:
      try:
        line = self._filehandle.next()
      except StopIteration:
        raise FormatError("Invalid VCF: No content")
      self._line_num+=1
      line = line.rstrip('\r\n')

      if line[0] == '#':
        if line[1] == '#':
          self._meta_header += line+'\n'
        elif line[0:6].upper() == '#CHROM':
          self._column_header += line+'\n'
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

  def get_meta_header(self):
    return self._meta_header

  def get_column_header(self):
    return self._column_header

  def get_header(self):
    return self._meta_header + self._column_header

  def get_sample_names(self):
    return self._sample_names

  def set_sample_names(self, sample_names):
    if isinstance(sample_names, list):
      self._sample_names = sample_names
      return True
    else:
      return False


#TODO: Separate systems for saying an attribute is not initialized and that its
#      value in the file is null
#TODO: Replace getters, setters with property decorators
class VCFSite(object):

  def __init__(self):
    # internal use
    self._initialized = False
    self._sample_names = None
    self._line_num = '[N/A]'
    self._modified = [True] * 10
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
    self._modified = [False] * 10

  # Attributes are unset until asked for
  def get_chrom(self):
    if self._chrom is None and self._initialized:
      if self._columns[0] == '.':
        raise FormatError("Invalid VCF: no CHROM in line "+self._line_num)
      self._chrom = self._columns[0]
    return self._chrom

  def get_pos(self):
    if self._pos is None and self._initialized:
      if self._columns[1] == '.':
        raise FormatError("Invalid VCF: no POS in line "+self._line_num)
      try:
        self._pos = int(self._columns[1])
      except ValueError:
        raise FormatError("Invalid VCF: non-integer POS in line "
          +self._line_num)
    return self._pos

  def get_id(self):
    if self._id is None and self._initialized:
      if self._columns[2] == '.':
        self._id = ''
      else:
        self._id = self._columns[2]
    return self._id

  def get_ref(self):
    if self._ref is None and self._initialized:
      if self._columns[3] == '.':
        raise FormatError("Invalid VCF: no REF in line "+self._line_num)
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
        self._qual = ''
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
        self._filter = []
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
      self._set_initialized()
      return True
    else:
      return False

  def set_pos(self, pos):
    if isinstance(pos, int):
      self._pos = pos
      self._modified[1] = True
      self._set_initialized()
      return True
    else:
      return False

  def set_id(self, id):
    if isinstance(id, basestring) or id is None:
      self._id = id
      self._modified[2] = True
      self._set_initialized()
      return True
    else:
      return False

  def set_ref(self, ref):
    if isinstance(ref, basestring) or ref is None:
      self._ref = ref
      self._modified[3] = True
      self._set_initialized()
      return True
    else:
      return False

  #TODO: Alter the NVC data in the INFO and genotypes to be consistent with new
  #      list of ALTs
  def set_alt(self, alt):
    if isinstance(alt, list) or alt is None:
      self._alt = alt
      self._modified[4] = True
      self._set_initialized()
      return True
    else:
      return False

  def set_qual(self, qual):
    if isinstance(qual, (int, long, float)) or qual == '':
      self._qual = qual
      self._modified[5] = True
      self._set_initialized()
      return True
    else:
      return False

  def set_filter(self, filter):
    if isinstance(filter, list) or filter is None or filter is True:
      self._filter = filter
      self._modified[6] = True
      self._set_initialized()
      return True
    else:
      return False

  def set_info(self, info):
    if isinstance(info, dict) or info is None:
      self._info = info
      self._modified[7] = True
      self._set_initialized()
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
      self._set_initialized()
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


  def alt_to_variant(self, alt):
    # N.B.: the reference allele is handled by the SNV logic
    ref = self.get_ref()
    diff = len(alt) - len(ref)
    if diff < 0:    # deletion
      return 'd'+str(-diff)
    elif diff > 0:  # insertion
      return alt[:diff+1]
    else:           # SNV
      return alt[0]
    #TODO: support for complex substitutions? (e.g. ref "GAT", alt "CCG")


  def split(self):
    """Split a multi-sample VCFSite into multiple VCFSites, one per sample."""

    sites = []
    genotypes = copy.deepcopy(self.get_genotypes())
    for sample_name in self._sample_names:
      site = VCFSite()
      site.set_chrom(self.get_chrom())
      site.set_pos(self.get_pos())
      site.set_id(self.get_id())
      site.set_ref(self.get_ref())
      site.set_alt(copy.deepcopy(self.get_alt()))
      site.set_qual(self.get_qual())
      site.set_filter(copy.deepcopy(self.get_filter()))
      site.set_info(copy.deepcopy(self.get_info()))
      site.set_genotypes({sample_name:genotypes[sample_name]})
      sites.append(site)

    return sites


  def __str__(self):
    """Returns a VCF line with the site's current data (no newline)"""
    if True not in self._modified:
      return '\t'.join(self._columns)

    line = []
    line.append(self.get_chrom())
    line.append(self.get_pos())
    line.append(self.get_id() if self.get_id() else '.')
    line.append(self.get_ref())
    line.append(','.join(self.get_alt()) if self.get_alt() else '.')
    line.append(self.get_qual() if self.get_qual() else '.')
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


  def _set_initialized(self):
    has_line = bool(
      self._columns is not None and
      isinstance(self._columns, list) and
      len(self._columns) >= 10
    )
    fields_full = bool(
      self._sample_names is not None and
      self._chrom is not None and
      self._pos is not None and
      self._id is not None and
      self._ref is not None and
      self._alt is not None and
      self._qual is not None and
      self._filter is not None and
      self._info is not None and
      self._genotypes is not None
    )
    self._initialized = has_line or fields_full
    return self._initialized


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
