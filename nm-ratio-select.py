#! /usr/bin/python
# original author: Boris Rebolledo-Jaramillo
# maintainer: Nick Stoler
from __future__ import division
import sys
import pysam
import argparse

DESCRIPTION = "Filter a BAM by NM-tag"
OPT_DEFAULTS = {'threshold':2.0}
parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.set_defaults(**OPT_DEFAULTS)
parser.add_argument('input', metavar='alignment.bam',
  help='Input filehame (BAM file). N.B.: Must be sorted by read name!')
parser.add_argument('-o', '--output', metavar='new-align.bam',
  help='Output filehame. Overrides default of prepending input filename with '
    '"nm-ratio.".')
parser.add_argument('-t', '--threshold', metavar='pct', type=float,
  help='NM edit distance threshold. In percentage per bp of read length.')
args = parser.parse_args()

if args.output:
  output = args.output
else:
  output = "nm-ratio."+args.input

def get_nm(read):
  tags = dict(read.tags)
  return float(tags['NM'])

sam = pysam.Samfile(args.input, 'rb')
out = pysam.Samfile(output, 'wb', template=sam)

for read in sam:
  try:
    read1 = read
    read2 = sam.next()
    if (get_nm(read1) <= (read1.rlen * args.threshold/100) and
        get_nm(read2) <= (read2.rlen * args.threshold/100)):
      out.write(read1)
      out.write(read2)
    else:
      pass

  except StopIteration:
    sam.close()
    out.close()
