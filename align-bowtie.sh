#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

DefaultChunkMbs=512
RefdirDefault=refdir
RequiredCommands='bowtie bowtie-build samtools awk'

Usage="Usage: \$ $(basename $0) [options] [-b align.bam | -s align.sam] ref.fa reads_1.fq reads_2.fq
-t: Number of threads for bowtie and bowtie-build to use (default: 1).
-m: Number to pass to bowtie's --chunkmbs option (default: $DefaultChunkMbs).
-c: Do extra cleanup afterward: Delete any reference index files created by this script."

function main {

  # Read in arguments and check them.
  sam=
  bam=
  cleanup=
  threads=1
  chunkmbs="$DefaultChunkMbs"
  while getopts "hs:b:m:t:c" opt; do
    case "$opt" in
      s) sam="$OPTARG";;
      b) bam="$OPTARG";;
      t) threads="$OPTARG";;
      m) chunkmbs="$OPTARG";;
      c) cleanup=true;;
      [h?]) fail "$Usage";;
    esac
  done
  # Get positional arguments.
  ref=${@:$OPTIND:1}
  fq1=${@:$OPTIND+1:1}
  fq2=${@:$OPTIND+2:1}

  # Validate arguments.
  if ! [[ "$bam" ]] && ! [[ "$sam" ]]; then
    fail "$Usage"$'\n'$'\n'"Error: Must provide an output sam or bam name."
  elif [[ "$bam" ]] && [[ -e "$bam" ]]; then
    fail "Error: Output bam \"$bam\" already exists."
  elif [[ "$sam" ]] && [[ -e "$sam" ]]; then
    fail "Error: Output sam \"$sam\" already exists."
  fi
  for path in "$ref" "$fq1" "$fq2"; do
    if ! [[ -s "$path" ]]; then
      fail "Error: Input file nonexistent or empty: $path"
    fi
  done

  # Determine names of the SAM and BAM files.
  if [[ "$sam" ]]; then
    format=sam
  elif [[ "$bam" ]]; then
    format=bam
    sam="$(echo "$bam" | sed -E 's/\.[^.]+$//').sam"
    if [[ -e "$sam" ]]; then
      fail "Error: Intermediate file \"$sam\" already exists."
    fi
  fi

  # Check for required commands.
  for cmd in $RequiredCommands; do
    if ! which "$cmd" >/dev/null 2>/dev/null; then
      fail "Error: command \"$cmd\" not found."
    fi
  done

  # Check version of bowtie-build.
  # Only version 1.2.1 and above had --threads option.
  indexer_is_threaded=$(bowtie-build --version | awk '
    $1 == "bowtie-build" && $2 == "version" {
      split($3, fields, ".")
      maj_min = fields[1] "." fields[2]
      if (maj_min > 1.2) {
        print "yes"
      } else if (maj_min == 1.2 && fields[3] >= 1) {
        print "yes"
      }
    }')
  if [[ $indexer_is_threaded ]]; then
    indexer_threads="--threads $threads"
  else
    indexer_threads=
  fi

  echo "\
ref:    $ref
fastq1: $fq1
fastq2: $fq2
sam:    $sam
bam:    $bam" >&2

  # Does the reference look indexed?
  index_missing=
  declare -a to_be_cleaned
  i=0
  for ext in 1.ebwt 2.ebwt 3.ebwt 4.ebwt rev.1.ebwt rev.2.ebwt; do
    if ! [[ -s "$ref.$ext" ]]; then
      index_missing=true
      if [[ "$cleanup" ]]; then
        to_be_cleaned["$i"]="$ref.$ext"
        i=$((i+1))
      fi
    fi
  done

  # Perform alignment.
  if [[ "$index_missing" ]]; then
    bowtie-build -f $indexer_threads --offrate 1 "$ref" "$ref" >/dev/null
  fi
  bowtie --chunkmbs "$chunkmbs" --threads "$threads" -f --sam -v 3 \
    "$ref" -1 "$fq1" -2 "$fq2" "$sam"

  # Check output and convert if requested.
  if [[ "$format" == sam ]]; then
    if [[ -s "$sam" ]]; then
      echo "Success. Output located in \"$sam\"." >&2
    else
      fail "Warning: Output file \"$bam\" empty or missing."
    fi
  elif [[ "$format" == bam ]]; then
    samtools view -Sb "$sam" | samtools sort -o - dummy > "$bam"
    if [[ -s "$bam" ]]; then
      samtools index "$bam"
      to_be_cleaned["${#to_be_cleaned[@]}"]="$sam"
      echo "Success. Output located in \"$bam\"." >&2
    else
      fail "Warning: Output file \"$bam\" empty or missing."
    fi
  fi

  # Cleanup.
  i=0
  while [[ "$i" -lt "${#to_be_cleaned[@]}" ]]; do
    if [[ -f "${to_be_cleaned[$i]}" ]]; then
      rm "${to_be_cleaned[$i]}"
    fi
    i=$((i+1))
  done
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
