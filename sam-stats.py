#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import collections
import samreader
import samflags

OPT_DEFAULTS = {'str':'string', 'int':0}
USAGE = """$ %(prog)s [options] alignments.sam
       $ samtools view -h alignments.bam | %(prog)s [options]"""
DESCRIPTION = """Print stats about the alignments in a SAM file."""
EPILOG = """"""

def main():

  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('file', metavar='alignments.sam', nargs='?',
    help='Input SAM file.')
  parser.add_argument('-s', '--name-sorted', action='store_true',
    help='The input SAM is sorted by read name. Allows reporting information '
      'on the paired status of the reads (assumes paired reads have the same '
      'name).')

  args = parser.parse_args()

  if args.file:
    filehandle = open(args.file)
  else:
    filehandle = sys.stdin

  flag_totals = collections.OrderedDict()
  for flag in samflags.FLAGS:
    flag_totals[flag] = 0

  last_name = None
  reads_set = []
  pair_stats = {'pairs':0, 'singletons':0, 'extras':0}
  for read in samreader.read(filehandle):
    # Add flags to totals
    read['flags'] = samflags.decompose(int(read['flag']))
    for flag in flag_totals:
      if read['flags'].get(flag):
        flag_totals[flag]+=1
    # Add pair stats, if sorted by name
    if args.name_sorted:
      name = read['qname']
      if last_name is not None and name != last_name:
        # We've started a new set of reads
        pair_stats = read_pair_stats(reads_set, pair_stats)
        reads_set = []
      reads_set.append(read)
      last_name = name

  print "flags totals:"
  for (flag, total) in flag_totals.items():
    width = max(len(str(flag)), len(str(total)))
    sys.stdout.write((" {:"+str(width)+"}").format(flag))
  print
  for (flag, total) in flag_totals.items():
    width = max(len(str(flag)), len(str(total)))
    sys.stdout.write((" {:"+str(width)+"}").format(total))
  print "\n"

  if args.name_sorted:
    pair_stats = read_pair_stats(reads_set, pair_stats)
    for (category, total) in pair_stats.items():
      print "{:11} {}".format(category+":", total)

  if filehandle is not sys.stdin:
    filehandle.close()


def read_pair_stats(reads_set, pair_stats):
  (paired_reads, single_reads, failed_reads) = partition_reads(reads_set)
  num_paired = len(paired_reads)
  num_single = len(single_reads)
  num_failed = len(failed_reads)
  # found two paired reads?
  if num_paired == 2:
    pair_stats['pairs'] += 1
  elif num_paired != 0:
    raise Exception('Uneven of "paired" reads: {} "paired" reads named '
      '"{}".'.format(num_paired, paired_reads[0]['qname']))
  # found one singleton read?
  if num_single == 1:
    pair_stats['singletons'] += 1
  elif num_single != 0:
    raise Exception('>1 "singleton" reads with the same name: {} "singleton" '
      'reads named "{}".'.format(num_single, single_reads[0]['qname']))
  # should only find paired reads OR singleton reads, not both
  if num_paired > 0 and num_single > 0:
    raise Exception('Found both "paired" reads and "singleton" reads for the '
      'same name. {} "paired" reads and {} "singleton" reads named '
      '"{}".'.format(num_paired, num_single, num_paired[0]['qname']))
  pair_stats['extras'] += num_failed
  return pair_stats


def partition_reads(reads_set):
  """Split a list of reads into pairs, singletons, and ones that failed filters.
  Failed reads are unmapped, secondary alignments, supplementary alignments, or
  "mapped in proper pair" but don't give a valid location for its mate.
  Paired reads are ones that didn't fail and are properly paired.
  Single reads are ones that didn't fail and aren't properly paired.
  Returns 3 lists: (paired_reads, single_reads, failed_reads)"""
  paired_reads = []
  single_reads = []
  failed_reads = []
  # Which reads pass the filters?
  for read in reads_set:
    flags = read['flags']
    # is unmapped or secondary/supplementary alignment?
    if flags[4] or flags[256] or flags.get(2048):
      failed_reads.append(read)
    # doesn't give location of mate?
    elif (read['rnext'] == '*' or read['pnext'] == '0' or
        read['rnext'] == ''  or read['pnext'] == ''):
      failed_reads.append(read)
    else:
      paired_reads.append(read)
  # Double-check if passing reads are properly paired or singletons
  if len(paired_reads) == 1:
    single_reads.append(paired_reads[0])
    del(paired_reads[0])
  else:
    # Move reads that aren't properly paired to the "failed" list
    i = 0
    while i < len(paired_reads):
      if not flags[2]:
        failed_reads.append(paired_reads[i])
        del(paired_reads[i])
      else:
        i+=1
  return (paired_reads, single_reads, failed_reads)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()
