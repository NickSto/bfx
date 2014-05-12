#!/usr/bin/env python
# original author: Nick Stoler
import os
import sys
import argparse

USAGE = "$ %(prog)s [options] flagint [flagint [flagint [..]]]"
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

  parser = argparse.ArgumentParser(
    usage=USAGE, description=DESCRIPTION, epilog=EPILOG
  )
  parser.add_argument('flags', metavar='flagint', nargs='+',
    help='The integer form of the flags for the read (e.g. "83" for 1, 2, 16, '
      'and 64).')
  parser.add_argument('-s', '--only-set', action='store_true',
    help='Only print the set flags.')

  args = parser.parse_args()

  for flag in args.flags:
    if len(args.flags) > 1:
      print flag+":"
    try:
      flagint = int(flag)
    except ValueError:
      fail('Error: flag must be an integer. Failed on "'+flag+'"')
    if flagint > 2**NUM_FLAGS - 1:
      fail('Error: invalid flag (greater than '+str(2**NUM_FLAGS - 1)
        +'. Failed on "'+flag+'"')

    print_flags(get_set_flags(flagint), print_unset=not(args.only_set))


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
        print "[X] %4d %s" % (2**i, MEANINGS[2**i].upper())
      else:
        print "[ ] %4d %s" % (2**i, MEANINGS[2**i].lower())
    else:
      if set_flags[i]:
        print "%4d %s" % (2**i, MEANINGS[2**i])


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
