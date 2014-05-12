#!/usr/bin/env bash
# original author: Nick Stoler
set -ue

command=$(basename $0)
USAGE="Usage: \$ $command region seq.fa [seq2.fa [seq3.fa [...]]]
       \$ cat seq.fa | $command region
  e.g. \$ $command chrM:1020-1080 hg19.fa
Use bioawk to print the sequence at the given coordinates in a FASTA (or FASTQ)
file. Use UCSC format for \"region\". Chromosome and ending coordinate are
optional. If no ending coordinate, it prints a single base at the start
coordinate. If no chromosome, it prints for every chromosome in the file.
WARNING:
bioawk seems to load each chromosome entirely into memory. For hg19.fa, it uses
about 1GB."

if [[ $# -lt 1 ]] || [[ $1 == '-h' ]]; then
  command=$(basename $0)
  echo "$USAGE" >&2
  exit 1
fi
region=$1

if ! which bioawk >/dev/null 2>/dev/null; then
  echo "Error: cannot find bioawk on PATH" >&2
  exit 1
fi

if expr index $region ':' >/dev/null; then
  colon=$(expr index $region ':')
  chrom=${region:0:$((colon-1))}
  region=${region:colon}
  restrict='$name == "'$chrom'" '
else
  restrict=''
fi

if expr index $region '-' >/dev/null; then
  dash=$(expr index $region '-')
  start=${region:0:$((dash-1))}
  end=${region:dash}
  length=$((end-start+1))
else
  start=$region
  length=1
fi

script=$restrict'{print substr($seq, '$start', '$length')}'
if [[ $# -gt 1 ]]; then
  for file in ${@:2}; do
    bioawk -c fastx "$script" $file
  done
else
  bioawk -c fastx "$script"
fi
