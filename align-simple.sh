#!/usr/bin/env bash
# original author: Nick Stoler
if [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue
PLATFORM=${PLATFORM:=ILLUMINA}
BWA_OPTS=${BWA_OPTS:=-M -t 16}
if [[ $# -le 1 ]]; then
  echo "Usage: $(basename $0) ref.fa reads_1.fq reads_2.fq out.bam [sample_id]"
  exit 1
fi
# Read arguments
read ref fastq1 fastq2 bam sample <<< "$@"
base=$(echo "$bam" | sed 's/\.bam$//')
if ! [[ $sample ]]; then
  sample=$(basename $base)
fi
# Index reference
if ! [[ -s $ref.amb ]] || ! [[ -s $ref.ann ]] || ! [[ -s $ref.bwt ]] || \
    ! [[ -s $ref.sa ]] || ! [[ -s $ref.pac ]]; then
  algo="bwtsw"
  if [[ $(du -sb $ref | awk '{print $1}') -lt 2000000000 ]]; then
    algo="is"
  fi
  bwa index -a $algo $ref
fi
# Align, adding read group
bwa mem $BWA_OPTS -R "@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM" $ref $fastq1 $fastq2 > $base.sam
samtools view -Sb $base.sam > $base.tmp.bam
samtools sort $base.tmp.bam $base
samtools index $bam
rm $base.sam $base.tmp.bam
