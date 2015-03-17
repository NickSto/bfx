#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

function main {

  # Code from Pierre Lindenbaum: https://www.biostars.org/p/6544/#6562
  paste f1.fastq f2.fastq |\ #merge the two fastqs
    awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\ #merge by group of 4 lines
    shuf  |\ #shuffle
    head |\ #only 10 records
    sed 's/\t\t/\n/g' |\ #restore the delimiters
    awk '{print $1 > "file1.fastq"; print $2 > "file2.fatsq"}' #split in two files.

}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
