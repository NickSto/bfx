#!/usr/bin/env python
# original author: Nick Stoler
import sys
import argparse
import samflags

USAGE = "$ %(prog)s [options] flagint [flagint [flagint [..]]]"
DESCRIPTION = """Decompose a SAM flag integer and print the individual flags that are set."""
EPILOG = """Meanings taken from http://picard.sourceforge.net/explain-flags.html but double-checked
against the SAM spec."""

MAX_VALUE = 2**len(samflags.FLAGS) - 1

def main():

  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)
  parser.add_argument('flags', metavar='flagint', nargs='+',
    help='The integer form of the flags for the read (e.g. "83" for 1, 2, 16, and 64).')
  parser.add_argument('-s', '--only-set', action='store_true',
    help='Only print the set flags.')

  args = parser.parse_args()

  for flag in args.flags:
    if len(args.flags) > 1:
      print flag+":"
    try:
      flagint = int(flag)
    except ValueError:
      fail('Error: flag must be an integer. Failed on "%s".' % flag)
    if flagint > MAX_VALUE:
      fail('Error: invalid flag (greater than %d). Failed on "%s".' % (MAX_VALUE, flag))

    print_flags(samflags.decompose(flagint), print_unset=not(args.only_set))


def print_flags(flags_set, print_unset=False):
  for (flag, meaning) in samflags.FLAGS.items():
    if print_unset:
      if flags_set[flag]:
        print "[X] %4d %s" % (flag, meaning.upper())
      else:
        print "[ ] %4d %s" % (flag, meaning.lower())
    else:
      if flags_set[flag]:
        print "%4d %s" % (flag, meaning)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()
