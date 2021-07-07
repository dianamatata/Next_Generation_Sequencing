#!/bin/bash
# forward read = SRR519926_1.fastq , reverse read SRR519926_2.fastq
# how do I check the supposed amount of reads? https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR519926
# 400,596 spots in the page
# echo $(cat SRR519926_1.fastq | wc -l)/4 | bc
# 1 read = 4lines

for file in reads/*.fastq; do
	reads=$(echo $(cat $file | wc -l)/4 | bc)
	echo "$file : $reads reads"
done
