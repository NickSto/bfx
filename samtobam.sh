#!/bin/bash
# original author: Nick Stoler
set -ue

if [ ${#@} -lt 1 ]; then
  echo "Produces a sorted, indexed BAM file from a SAM file." 
  echo "USAGE:"
  echo "  $ samtobam.sh align.sam"
  echo -n "Note: it will try to get the basename of the file by removing any "
  echo "ending '.sam'"
  echo "So the above command will produce 'align.bam'"
  exit 0
fi

sam=$1
base=$(echo "$sam" | sed -r 's/\.sam$//')
tmpbam=$base.tmp.bam
bam=$base.bam

set -x

samtools view -Sb $sam > $tmpbam

samtools sort $tmpbam $base

samtools index $bam

rm $tmpbam
