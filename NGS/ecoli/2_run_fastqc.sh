#!/bin/bash
# alternative: fastqc *.fastq
# fastqc creates html file with results and .zip files

for file in reads/*.fastq; do
	fastqc $file
done
