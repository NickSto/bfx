#!/usr/bin/env python
from __future__ import print_function
import sys
import os

script_name = os.path.basename(sys.argv[0])
quiet = False

USAGE = script_name+""" reads.fq

This helps determine what quality score encoding is being used in a FASTQ file.
It will print the highest and lowest quality character found.
Warning: this does not work with FASTQ files containing multi-line sequences."""

def main():

  if len(sys.argv) > 1:
    filepath = sys.argv[1]
  else:
    fail("Please provide input FASTQ filename.")

  line_num = 0
  read_num = 0
  ord_min = 128
  ord_max = 0
  seqlines = 0
  with open(filepath, 'rU') as filehandle:
    for line in filehandle:
      line = line.strip()
      line_num+=1
      line_type = line_num % 4
      # 1st ID line
      if line_type == 1:
        read_num+=1
        if line[0] != '@':
          fail("Format error on line "+str(line_num))
      # sequence line
      elif line_type == 2:
        pass
      # 2nd ID line
      elif line_type == 3:
        if line[0] != '+':
          fail("Format error on line "+str(line_num))
      # quality line
      elif line_type == 0:
        if quiet:
          (ord_min, ord_max) = minmax(line, ord_min, ord_max)
        else:
          (ord_min, ord_max) = minmax(line, ord_min, ord_max, read_num)

  print("Lowest:  "+str(ord_min)+" ("+chr(ord_min)+")")
  print("Highest: "+str(ord_max)+" ("+chr(ord_max)+")")


def minmax(line, ord_min, ord_max, read_num=None):
  for char in line:
    if ord(char) < ord_min:
      ord_min = ord(char)
      if read_num is not None:
        print("smaller:  "+char+"         read "+str(read_num))
    if ord(char) > ord_max:
      ord_max = ord(char)
      if read_num is not None:
        print("bigger:        "+char+"    read "+str(read_num))
  return (ord_min, ord_max)


def fail(message):
  sys.stderr.write(USAGE+"\n")
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  main()