# With samtools view you can easily filter your alignment file based on flags
# info: http://www.htslib.org/doc/samtools-view.html

cd ~/ecoli/alignment_output

# Filtering against unmapped reads (leaving only mapped reads) with samtools view would look like this:
# -F Do not output alignments with any bits set in INT present in the FLAG field. 
# -f Only output alignments with all bits set in INT present in the FLAG field.
# 0x4 flag for unmapped reads
samtools view -bh -F 0x4 SRR519926.sorted.bam > SRR519926.sorted.mapped.bam

#  How many reads are in there? -c for count
samtools view -f 0x4 -c SRR519926.sorted.bam

# decompose
samtools view -bh -f 0x4 SRR519926.sorted.bam > SRR519926.sorted.unmapped.bam
samtools view -c SRR519926.sorted.unmapped.bam

# is this consistent with the flagstat output?
total_reads=$(cat samtools_flagstat_output.txt | grep 'total' | cut -d' ' -f1)
reads_mapped=$(cat samtools_flagstat_output.txt | grep 'mapped' | cut -d' ' -f1 | head -1)
reads_unmapped=$(samtools view -f 0x4 -c SRR519926.sorted.bam)

# Filter our sorted and indexed BAM file for the region between 2000 and 2500 kb, and output it as a BAM file with a header.
cd ~/ecoli/ref_genome
grep ">" ecoli-strK12-MG1655.fasta
# >U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome
cd ~/ecoli/alignment_output
samtools view -bh SRR519926.sorted.bam U00096.3:2000000-2500000 > SRR519926.sorted.region.bam

