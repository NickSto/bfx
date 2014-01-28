#!/usr/bin/env python
# original author: Nick Stoler
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'only_set':False}
USAGE = "USAGE: %prog [options] flagint [flagint [flagint [..]]]"
DESCRIPTION = """Decompose a SAM flag integer and print the individual flags
that are set."""
EPILOG = """Meanings taken from http://picard.sourceforge.net/explain-flags.html
but double-checked against the SAM spec."""

MEANINGS = {
  1:'read paired',
  2:'read mapped in proper pair',
  4:'read unmapped',
  8:'mate unmapped',
  16:'read reverse strand',
  32:'mate reverse strand',
  64:'first in pair',
  128:'second in pair',
  256:'not primary alignment',
  512:'read fails platform/vendor quality checks',
  1024:'read is PCR or optical duplicate',
  2048:'supplementary alignment',
}
NUM_FLAGS = len(MEANINGS)

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  # parser.add_option('-s', '--str', dest='str',
  #   default=OPT_DEFAULTS.get('str'), help='default: %default')
  # parser.add_option('-i', '--int', dest='int', type='int',
  #   default=OPT_DEFAULTS.get('int'), help='')
  # parser.add_option('-f', '--float', dest='float', type='float',
  #   default=OPT_DEFAULTS.get('float'), help='')
  parser.add_option('-s', '--only-set', dest='only_set',
    action='store_const', const=not OPT_DEFAULTS.get('only_set'),
    default=OPT_DEFAULTS.get('only_set'),
    help='Print the unset flags too')

  (options, arguments) = parser.parse_args()

  if not arguments:
    parser.print_help()
    fail("Please provide a SAM flag.")

  for flag in arguments:
    if len(arguments) > 1:
      print flag+":"
    try:
      flagint = int(flag)
    except ValueError:
      fail('Error: flag must be an integer. Failed on "'+flag+'"')
    if flagint > 2**NUM_FLAGS - 1:
      fail('Error: invalid flag (greater than '+str(2**NUM_FLAGS - 1)
        +'. Failed on "'+flag+'"')

    print_flags(get_set_flags(flagint), print_unset=not(options.only_set))


def get_set_flags(flagint):
  flagbin = ('{0:0'+str(NUM_FLAGS)+'b}').format(flagint)
  i = 0
  set_flags = [False] * NUM_FLAGS
  for bit in reversed(flagbin):
    set_flags[i] = bit == '1'
    i+=1
  return set_flags


def print_flags(set_flags, print_unset=False):
  for i in range(len(set_flags)):
    if print_unset:
      if set_flags[i]:
        print "SET: %4d %s" % (2**i, MEANINGS[2**i].upper())
      else:
        print "not: %4d %s" % (2**i, MEANINGS[2**i].lower())
    else:
      if set_flags[i]:
        print "%4d %s" % (2**i, MEANINGS[2**i])


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
