#!/usr/bin/env bash
#Check out the statistics of the E. coli alignment by using samtools flagstat
samtools flagstat alignment_output/SRR519926.sam > alignment_output/samtools_flagstat_output.txt
echo -e '\nprint flagstat output\n'
cat alignment_output/samtools_flagstat_output.txt

# Of the reads, 38.44% is properly paired. The rest isn’t. Proper pairing is quite hard to interpret. It usually means that the 0x2 flag (each segment properly aligned according to the aligner) is false. In this case it means that the insert size is high for a lot of sequences. That is because the insert size distribution is very wide. You can find info on insert size distribution like this:

samtools stats alignment_output/SRR519926.sam | grep ^SN | cut -f 2,3 > alignment_output/stats.txt
echo -e '\nprint stats\n'
cat alignment_output/stats.txt

# Now look at insert size average and insert size standard deviation
cat alignment_output/stats.txt | grep 'insert size'

# The command samtools view is very versatile. It takes an alignment file and writes a filtered or processed alignment to the output. You can for example use it to compress your SAM file into a BAM file.
# -h for header -b for output bam
samtools view -bh alignment_output/SRR519926.sam > alignment_output/SRR519926.bam

# To look up specific alignments, it is convenient to have your alignment file indexed. An indexing can be compared to a kind of ‘phonebook’ of your sequence alignment file. Indexing can be done with samtools as well, but it first needs to be sorted on coordinate (i.e. the alignment location)
samtools sort alignment_output/SRR519926.bam > alignment_output/SRR519926.sorted.bam
samtools index alignment_output/SRR519926.sorted.bam

#  How much was the disk space reduced by compressing the file?
cd alignment_output
ls -lh * | grep 'SRR519926.sam'
ls -lh * | grep 'SRR519926.bam'
