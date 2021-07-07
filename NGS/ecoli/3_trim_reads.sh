##!/usr/bin/env bash

# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

TRIMMED_DIR=~/ecoli/trimmed_data
READS_DIR=~/ecoli/reads
ADAPTERS=/opt/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/TruSeq3-PE.fa

mkdir $TRIMMED_DIR

trimmomatic \
PE \
-threads 1 \
-phred33 \
$READS_DIR/SRR519926_1.fastq \
$READS_DIR/SRR519926_2.fastq \
$TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
$TRIMMED_DIR/unpaired_trimmed_SRR519926_1.fastq \
$TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
$TRIMMED_DIR/unpaired_trimmed_SRR519926_2.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 \
SLIDINGWINDOW:4:5 \
LEADING:5 \
TRAILING:5 \
MINLEN:25
