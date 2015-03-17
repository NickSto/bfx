#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

PCT_DEFAULT=1

USAGE="Usage: \$ $(basename $0) [options] input_1.fq input_2.fq output_1.fq output_2.fq
Options:
-p: Proportion of input reads to output, in percent. Default: $PCT_DEFAULT%.
-n: Absolute number of read pairs to output.
WARNING: This assumes one-line reads (each record is 4 lines)."

function main {

  # read in arguments
  pct="$PCT_DEFAULT"
  pairs=""
  while getopts ":p:n:h" opt; do
    case "$opt" in
      p) pct="$OPTARG";;
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

  # Check lengths of fastq files.
  wc1=$(wc -l $in1 | awk '{print $1}')
  wc2=$(wc -l $in2 | awk '{print $1}')
  if [[ $wc1 == $wc2 ]]; then
    total_pairs=$((wc1/4))
    if [[ $((total_pairs*4)) != $wc1 ]]; then
      fail "Error: Input fastq files \"$in1\" and \"$in2\"
have line counts that aren't multiples of 4. Make sure they don't contain multi-
line sequences."
    fi
  else
    fail "Error: Input fastq files \"$in1\" and \"$in2\"
are different lengths: $wc1 and $wc2 lines, respectively."
  fi

  # Convert percent of reads to absolute number of reads.
  if ! [[ $pairs ]]; then
    pairs=$(printf %.0f $(echo "$total_pairs*$pct/100" | bc -l))
  fi

  # Code from Pierre Lindenbaum: https://www.biostars.org/p/6544/#6562
  # Merge the two fastqs.
  paste $in1 $in2 |\
    # Merge by groups of 4 lines, print all 8 (for both pairs) on one line
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
    # Shuffle lines.
    shuf  |\
    # Cut off after $pairs lines.
    head -n $pairs |\
    # Restore the delimiters.
    sed 's/\t\t/\n/g' |\
    # Split in two files.
    awk -F '\t' '{print $1 > "'$out1'"; print $2 > "'$out2'"}'

}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
