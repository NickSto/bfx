#!/usr/bin/env bash
# original author: Nick Stoler
# Automates all the steps to generate an indexed BAM file using BWA and SAMTools
set -ue

REQUIRED="bwa samtools" # must be in PATH
REF_THRES=2000000000
BWA_OPTS_DEFAULT="-M"

USAGE="USAGE:
  \$ $(basename $0) [-d outdir] [-o \"bwa opts\"] ref.fa reads_1.fq [reads_2.fq]
Align FASTQ with BWA-MEM and generate an indexed BAM file.
This names the output BAM using the base derived from the first fastq file,
minus any .fq or .fastq extension, and minus any _1 or _2 suffix. It indexes the
reference if it hasn't been already, choosing bwtsw unless the file is smaller
than $REF_THRES bytes. It will not overwrite any existing SAM or BAM file.
Options:
-d outdir:     The output directory. The BAM and all other new files (SAM, BAI)
               will be created here. Default is the directory containing the
               first FASTQ file.
-o \"bwa opts\": The options you want to pass to the \"bwa mem\" command. Make
               sure to enclose your options in quotes.
               Default: \"$BWA_OPTS_DEFAULT\""

function fail {
  echo "$@"
  exit 1
}
function echev {
  echo "+ $1"
  eval "$1"
}

# Check that it can find tools it needs like bwa and samtools
for cmd in $REQUIRED; do
  if ! which $cmd >/dev/null; then
    fail "Error: could not find $cmd in PATH"
  fi
done

########## Gather settings ##########

# Get options
outdir=""
bwa_opts="$BWA_OPTS_DEFAULT"
while getopts ":d:o:h" opt; do
  case "$opt" in
    d) outdir=$OPTARG;;
    o) bwa_opts="$OPTARG";;
    h) fail "$USAGE";;
  esac
done

# get positional arguments
narg=$OPTIND
while [[ $narg -le $# ]]; do
  arg=${@:$narg:1}
  if [[ ${arg:0:1} == '-' ]]; then
    fail "Error: options like $arg must come before positional arguments."
  fi
  narg=$((narg+1))
done
positionals=$((narg-OPTIND))
if [[ $positionals -lt 2 ]]; then
  fail "$USAGE"
fi
ref="${@:$OPTIND:1}"
fq1="${@:$OPTIND+1:1}"
fq2="${@:$OPTIND+2:1}"

# set $outdir, check it exists
if [[ -z $outdir ]]; then
  outdir=$(dirname $fq1)
fi
if [[ ! -d $outdir ]]; then
  fail "Error: could not find output directory $outdir"
fi

# get basename from fastq filename
# remove .fq
basename=$(basename $fq1 .fq)
# if that doesn't work, try .fastq
if [[ $basename == $fq1 ]]; then
  basename=$(basename $fq1 .fastq)
fi
# remove _1 or _2
if [[ $basename =~ _[12]$ ]]; then
  lastchar=$((${#basename}-2))
  basename=${basename:0:$lastchar}
fi

# create rest of filenames, print settings
bamtmp=$outdir/$basename.tmp.bam
bam=$outdir/$basename.bam

echo "\
Settings:
  basename: $basename
  ref:      $ref
  fq1:      $fq1"
if [[ -n $fq2 ]]; then
  echo "  fq2:      $fq2"
fi
echo "\
  bamtmp:   $bamtmp
  bam:      $bam
  bwa opts: $bwa_opts"

# check that input files exist
files="$ref $fq1 $fq2"
for file in $files; do
  if [[ ! -s $file ]]; then
    fail "Error: File $file nonexistent, inaccessible, or empty."
  fi
done
# stop if any of the output files exists already
for file in "$bam" "$bamtmp"; do
  if [[ -e "$file" ]]; then
    fail "Error: Output file $file already exists. Please rename either the \
existing file or the input fastq file."
  fi
done


########## Do alignment ##########

# does the reference look indexed?
if ! ( [[
      -s $ref.amb &&
      -s $ref.ann &&
      -s $ref.bwt &&
      -s $ref.pac &&
      -s $ref.sa
    ]] ); then
  # use bwtsw (safe for any size) unless it's definitely smaller than 2GB
  algo="bwtsw"
  if [[ $(du -sb $ref | awk '{print $1}') -lt $REF_THRES ]]; then
    algo="is"
  fi
  echev "bwa index -a $algo $ref"
fi

# alignment
echev "bwa mem $bwa_opts $ref $fq1 $fq2 | samtools view -Sb - > $bamtmp"

# sort BAM
echev "samtools sort -o $bamtmp dummy > $bam"

# index BAM
echev "samtools index $bam"

echev "rm $bamtmp"

echo -n "$fq1 "
if [[ -n $fq2 ]]; then
  echo -n "and $fq2 "
fi
echo "have been aligned to $ref into $bam"
