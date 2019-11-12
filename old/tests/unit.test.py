#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse

OPT_DEFAULTS = {}
USAGE = "%(prog)s [options] function.to.test"
DESCRIPTION = """Run test on a given function and print standardized results."""
EPILOG = """"""

REQUIREMENTS = {
  'LavReader':['lav'],
  'LavReader.convert':['lav'],
  'FastaLineGenerator.extract':['fasta','coord_file'],
}

def main():

  FUNCTIONS = {
    'LavReader':LavReader,
    'LavReader.convert':LavReader_convert,
    'FastaLineGenerator.extract':FastaLineGenerator_extract,
  }

  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('function', metavar='function.to.test',
    help="""The function to test. give one of: """+', '.join(FUNCTIONS))
  parser.add_argument('-l', '--lav',
    help="""Input LAV file, required for LavReader functions.""")
  parser.add_argument('-f', '--fasta',
    help="""Input FASTA file, required for FastaLineGenerator functions.""")
  parser.add_argument('-C', '--coord-file',
    help="""File containing input coordinates for FastaLineGenerator.extract
      (required). Coordinate format: chrom:start-end, chrom optional. E.g.
      "chr1:100-110", "300-400". One coordinate per line. Empty lines and
      comment lines (starting with #) are allowed.""")

  args = parser.parse_args()

  if args.function not in FUNCTIONS:
    fail('Error: function "%s" not supported. Please pick one from the list: %s'
      % (args.function, ', '.join(FUNCTIONS)))

  error = check_args(args, REQUIREMENTS)
  if error:
    parser.print_help()
    fail('\n'+error)

  FUNCTIONS[args.function](args)


##### UNIT TEST FUNCTIONS #####


def LavReader(args):
  import lavreader
  HIT_KEYS = ['filename', 'id', 'name', 'seqnum', 'revcomp', 'begin', 'end']
  SEGMENT_KEYS = ['begin', 'end']
  reader = lavreader.LavReader(args.lav, convert=False)
  for hit in reader:
    print '--- subject ---'
    for key in HIT_KEYS:
      print key+':\t'+str(hit.subject[key])
    print '--- query ---'
    for key in HIT_KEYS:
      print key+':\t'+str(hit.query[key])
    for alignment in hit:
      print '  score:',alignment.score
      print '  --- subject ---'
      for key in SEGMENT_KEYS:
        print '  '+key+':\t'+str(alignment.subject[key])
      print '  --- query ---'
      for key in SEGMENT_KEYS:
        print '  '+key+':\t'+str(alignment.query[key])
      for block in alignment:
        print '    identity:',block.identity
        print '    --- subject ---'
        for key in SEGMENT_KEYS:
          print '    '+key+':\t'+str(block.subject[key])
        print '    --- query ---'
        for key in SEGMENT_KEYS:
          print '    '+key+':\t'+str(block.query[key])


def LavReader_convert(args):
  import lavreader
  SEGMENT_KEYS = ['begin', 'end']
  reader = lavreader.LavReader(args.lav, convert=False)
  reader.convert()
  for hit in reader.hits:
    for alignment in hit.alignments:
      print '  --- subject ---'
      for key in SEGMENT_KEYS:
        print '  '+key+':\t'+str(alignment.subject[key])
      print '  --- query ---'
      for key in SEGMENT_KEYS:
        print '  '+key+':\t'+str(alignment.query[key])
      for block in alignment.blocks:
        print '    --- subject ---'
        for key in SEGMENT_KEYS:
          print '    '+key+':\t'+str(block.subject[key])
        print '    --- query ---'
        for key in SEGMENT_KEYS:
          print '    '+key+':\t'+str(block.query[key])


def FastaLineGenerator_extract(args):
  import fastareader
  coords = get_coords(args.coord_file)
  fasta = fastareader.FastaLineGenerator(args.fasta)
  for coord in coords:
    print fasta.extract(coord[1], coord[2], chrom=coord[0])


##### UTILITY FUNCTIONS #####

def get_coords(filepath):
  coords = []
  with open(filepath, 'rU') as filehandle:
    for raw_line in filehandle:
      line = raw_line.strip()
      if not line or line.startswith('#'):
        continue
      if ':' in line:
        (chrom, region) = line.split(':')
      else:
        chrom = None
        region = line
      assert '-' in region
      (start, end) = region.split('-')
      coords.append((chrom, int(start), int(end)))
  return coords


def check_args(args, all_requirements):
  """Returns an error message if a required option is missing, empty string
  otherwise.
  Includes a message about every missing requirement, not just the first."""
  errors = []
  requirements = all_requirements[args.function]
  for requirement in requirements:
    if getattr(args, requirement) is None:
      reqname = '--'+requirement.replace('_', '-')
      errors.append("Error: %s option required for function %s" %
        (reqname, args.function))
  return '\n'.join(errors)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
