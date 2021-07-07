#!/bin/bash
# alternative: fastqc *.fastq
# fastqc creates html file with results and .zip files

TRIMMED_DIR=~/ecoli/trimmed_data
for file in $TRIMMED_DIR/*.fastq; do
	fastqc $file
done
