#!/bin/bash
# Write a script that maps the reads with bowtie2 (see chapter 2 of read alignment), sorts them, takes only the mapped reads, and outputs them as a BAM file with a header.
# write in a pipe

ref=~/ecoli/ref_genome/ecoli-strK12-MG1655.fasta
TRIMMED_DIR=~/ecoli/trimmed_data
ALIGNED_DIR=~/ecoli/alignment_output

bowtie2 \
-x $ref \
-1 $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
-2 $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
| samtools sort - \
| samtools view -F 0x4 -bh - \
> $ALIGNED_DIR/SRR519926.sorted.mapped.frompipe.bam
# no indexing? it’s optional. However, for many downstream applications (like IGV), an index file (which you create with samtools index) is required.



