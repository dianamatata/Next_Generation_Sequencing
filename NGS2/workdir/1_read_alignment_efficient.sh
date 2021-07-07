#!/bin/bash

# piping and looping script.sh
# to use samtools in a pipe, the input argument needs to replaced with a -
# If the input file is absent, it reads from stdin. So, you could use it without a - replacing input or output files
# Also, some commands do not write by default to stdout, but to a specified file (this is the case for e.g. samtools fixmate and samtools markdup). In that case, also the output argument should be replaced with -.
# most commands, output -o files are optional. cmd writes in stdout

dir_data=/home/davalos/workdir/data
dir_align=/home/davalos/workdir/alignment_pipe
mkdir $dir_align

bwa index ${dir_data}/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
bwa index ${dir_data}/fastq/*.gz

for person in 'father' 'mother' 'son' ; do
        sample=${dir_align}/${person}

	bwa mem \
        ${dir_data}/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
        ${dir_data}/fastq/${person}_R1.fastq.gz \
        ${dir_data}/fastq/${person}_R2.fastq.gz \
	| samtools fixmate -m - - \
        | samtools sort \
        | samtools markdup - - \ 
	| samtools addreplacerg -r ID:${person} -r SM:${person} -r PL:ILLUMINA - \
	| samtools view -bh \
	> ${dir_align}/${person}.bam

	samtools index ${dir_align}/${person}.bam
done

# | samtools collate - \ By default the aligner will output the alignment collated, so you don’t need to run the command.
# markdup -s - - \ # s to print stats
