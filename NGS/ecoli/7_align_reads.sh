##!/usr/bin/env bash

TRIMMED_DIR=~/ecoli/trimmed_data
REFERENCE_DIR=~/ecoli/ref_genome/
ALIGNED_DIR=~/ecoli/alignment_output

mkdir $ALIGNED_DIR

bowtie2 \
-x $REFERENCE_DIR/ecoli-strK12-MG1655.fasta \
-1 $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
-2 $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
> $ALIGNED_DIR/SRR519926.sam
