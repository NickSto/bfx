#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

PAIRS_DEFAULT=10

USAGE="Usage: \$ $(basename $0) input_1.fq input_2.fq output_1.fq output_2.fq"

function main {

  # read in arguments
  pairs="$PAIRS_DEFAULT"
  while getopts ":n:h" opt; do
    case "$opt" in
      n) pairs="$OPTARG";;
      h) fail "$USAGE";;
    esac
  done
  in1="${@:$OPTIND:1}"
  in2="${@:$OPTIND+1:1}"
  out1="${@:$OPTIND+2:1}"
  out2="${@:$OPTIND+3:1}"
  if ! [[ "${@:$OPTIND+3:1}" ]]; then
    fail "$USAGE"
  fi

  if ! which shuf >/dev/null 2>/dev/null; then
    fail "Error: GNU \"shuf\" not found (BSD doesn't have it)."
  fi

  # Code from Pierre Lindenbaum: https://www.biostars.org/p/6544/#6562
  # merge the two fastqs
  paste $in1 $in2 |\
    # merge by group of 4 lines, print all eight (for both pairs) on one line
    awk '
      {
        printf("%s",$0);
        n++;
        if (n % 4 == 0) {
          printf("\n");
        } else {
          printf("\t\t");
        }
      }' |\
    # shuffle
    shuf  |\
    # only 10 records
    head -n $pairs |\
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
