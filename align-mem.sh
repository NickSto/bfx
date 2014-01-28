#!/usr/bin/env bash
# original author: Nick Stoler
# Automates all the steps to generate an indexed BAM file using BWA and SAMTools
# input is a reference genome file and a query reads file
set -ue

# must be in PATH
REQUIRED="bwa samtools"
REF_THRES=2000000000

if [[ $# -lt 2 ]]; then
  echo "USAGE:
  \$ $(basename $0) reference.fa reads_1.fq [reads_2.fq]
Align FASTQ with BWA-MEM and generate an indexed BAM file.
This names the output BAM using the base derived from the first fastq file,
minus any .fq or .fastq extension, and minus any _1 or _2 suffix. It indexes the
reference if it hasn't been already, choosing bwtsw unless the file is smaller
than $REF_THRES bytes. It will not overwrite any existing SAM or BAM file."
  exit 0
fi

for cmd in $REQUIRED; do
  if ! which $cmd >/dev/null; then
    echo "Error: could not find $cmd in PATH" 1>&2
    exit 1
  fi
done

basename=$(echo $2 | sed -E 's/\.f(ast)?q$//g')   # remove .fastq or .fq
basename=$(echo $basename | sed -E 's/_[12]$//g') # remove _1 or _2
ref=$1
fastq1=$2
fastq2=""
if [[ $# -gt 2 ]]; then
  fastq2=$3
fi
sam=$basename.sam
bamtmp=$basename.tmp.bam
bam=$basename.bam

echo "Using names:
  basename: $basename
  ref:      $ref
  fastq1:   $fastq1"
if [[ $fastq2 ]]; then
  echo "  fastq2:   $fastq2"
fi
echo "  sam:      $sam
  bamtmp:   $bamtmp
  bam:      $bam"

# stop if any of the permament output files exists already
existing=""
if [[ -e $sam ]]; then
  existing=$sam
elif [[ -e $bam ]]; then
  existing=$bam
fi
if [[ $existing ]]; then
  echo -n "Error: Output file $existing already exists. Please rename " 1>&2
  echo "either the existing file or the input fastq file." 1>&2
  exit 1
fi

# does the reference look indexed?
if ! ( [[
    -e $ref.amb &&
    -e $ref.ann &&
    -e $ref.bwt &&
    -e $ref.pac &&
    -e $ref.sa ]] ); then
  # use bwtsw (safe for any size) unless it's definitely smaller than 2GB
  algo="bwtsw"
  if [[ $(du -sb $ref | awk '{print $1}') -lt $REF_THRES ]]; then
    algo="is"
  fi
  cmd="bwa index -a $algo $ref"
  echo "\$ $cmd"; eval "$cmd"
fi

# alignment
cmd="bwa mem $ref $fastq1 $fastq2 > $sam"
echo "\$ $cmd"; eval "$cmd"

# make BAM
cmd="samtools view -Sb $sam > $bamtmp"
echo "\$ $cmd"; eval "$cmd"

# sort BAM
cmd="samtools sort $bamtmp $basename"
echo "\$ $cmd"; eval "$cmd"

# index BAM
cmd="samtools index $bam"
echo "\$ $cmd"; eval "$cmd"

cmd="rm $bamtmp"
echo "\$ $cmd"; eval "$cmd"

echo -n "$fastq1 "
if [[ $fastq2 ]]; then
  echo -n "and $fastq2 "
fi
echo "have been aligned to $ref into $bam"
