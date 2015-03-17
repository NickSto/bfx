#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

USAGE="Usage: \$ $(basename $0) input_1.fq input_2.fq output_1.fq output_2.fq"

function main {

  if [[ $# -lt 4 ]] || [[ $1 == '-h' ]]; then
    fail "$USAGE"
  fi

  in1="$1"
  in2="$2"
  out1="$3"
  out2="$4"

  # Code from Pierre Lindenbaum: https://www.biostars.org/p/6544/#6562
  # merge the two fastqs
  paste $in1 $in2 |\
    # merge by group of 4 lines
    awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
    # shuffle
    shuf  |\
    # only 10 records
    head |\
    # restore the delimiters
    sed 's/\t\t/\n/g' |\
    # split in two files.
    awk '{print $1 > "'$out1'"; print $2 > "'$out2'"}'

}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
