#!/usr/bin/env python
# original author: Boris Rebolledo-Jaramillo
# maintainer: Nick Stoler
from __future__ import division
import sys
import pysam
import argparse

DESCRIPTION = """Filter a BAM by NM-tag edit distance. Drop pairs of reads which
Have an edit distance greater than the threshold. Both reads must individually
pass the threshold. Unpaired reads are always dropped. Pairing is checked by
comparing read names."""
OPT_DEFAULTS = {'threshold':2.0}
parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.set_defaults(**OPT_DEFAULTS)
parser.add_argument('input', metavar='alignment.bam',
  help='Input BAM file. Must be sorted by read name!')
parser.add_argument('-o', '--output', metavar='new-align.bam',
  help='Output BAM file. Overrides default of prepending input filename with '
    '"nm-ratio.".')
parser.add_argument('-t', '--threshold', metavar='pct', type=float,
  help='NM edit distance threshold. In percentage per bp of read length.')
parser.add_argument('-F', '--unpaired-fatal', action='store_true',
  help='Fail on encountering an unpaired read.')
parser.add_argument('-P', '--no-pair-check', dest='check_pairs',
  action='store_false', default=True,
  help="Use original behaviour of not checking read names to make sure pairs "
    "are actual pairs. User must make sure no unpaired reads are present.")
args = parser.parse_args()

if args.output:
  output = args.output
else:
  output = "nm-ratio."+args.input

def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

def get_nm(read):
  tags = dict(read.tags)
  return float(tags['NM'])

sam = pysam.Samfile(args.input, 'rb')
out = pysam.Samfile(output, 'wb', template=sam)

total = 0
unpaired = 0
filtered = 0
for read1 in sam:
  total+=1
  try:
    read2 = sam.next()
    total+=1
    # check for actual pair
    if args.check_pairs:
      # step through reads, skipping unpaired ones
      while read1.qname != read2.qname:
        unpaired+=1
        if args.unpaired_fatal:
          fail("Error: unpaired read found:\n{}".format(read1.qname))
        read1 = read2
        read2 = sam.next()
        total+=1
        filtered+=1
    # check NM tag and print if the pair passes
    if (get_nm(read1) <= (read1.rlen * args.threshold/100) and
        get_nm(read2) <= (read2.rlen * args.threshold/100)):
      out.write(read1)
      out.write(read2)
    else:
      filtered+=2
  # StopIteration thrown on read2 = sam.next()
  except StopIteration:
    unpaired+=1
    if args.unpaired_fatal:
      fail("Error: unpaired read found:\n{}".format(read1.qname))
    sam.close()
    out.close()

# final report
if total == 0:
  fail("No reads found!")
elif args.check_pairs and unpaired/total > 0.5:
  print ('Very many "unpaired" reads found. Maybe the BAM wasn\'t sorted by '
    'read name?')
print 'total reads:    '+str(total)
if args.check_pairs:
  print 'unpaired reads: '+str(unpaired)
print 'reads removed:  '+str(filtered)
