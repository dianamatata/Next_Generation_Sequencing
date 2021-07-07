# ssh -i ~/Desktop/NGS2/key_davalos.pem davalos@18.158.121.17
#!/bin/bash
# conda activate variants

home_dir=/home/davalos/workdir
dir_data=$home_dir/data
dir_align=$home_dir/alignment/


# We’ll use bwa mem for the alignment. Like all alignment software, it requires an index of the reference genome. 
# http://bio-bwa.sourceforge.net/bwa.shtml

bwa index ${dir_data}/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
bwa index ${dir_data}/fastq/*.gz


# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.
# We’ll be using paired-end reads of three samples that can be found at ~/workdir/data/fastq. If we run bwa mem with default options, we need ref files and 2 paired reads

mkdir /home/davalos/workdir/alignment

for person in 'father' 'mother' 'son' ; do
	bwa mem \
	${dir_data}/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
	${dir_data}/fastq/${person}_R1.fastq.gz \
        ${dir_data}/fastq/${person}_R2.fastq.gz \
	> ${dir_align}/${person}.sam
done


# Check out the statistics of the alignment by using samtools flagstat.

for file in ~/workdir/alignment/*.sam ; do
	samtools flagstat $file > $file.stats
done


# compress your SAM file into a BAM file

for file in ${dir_align}/*.sam ; do
	name=$(echo $file | cut -d'.' -f1)
	samtools view -bh $name.sam > $name.bam
done


# check size difference with compression

ls -lh


# For variant analysis, it’s important to mark reads that possibly originated from PCR duplication. We can do that with samtools markdup. However, we can not directly run that on our .sam file nor on our compressed .bam file.
# samtools markdup – mark duplicate alignments in a coordinate sorted file

for person in 'father' 'mother' 'son' ; do
	
	name=${dir_align}/${person}
	# This first collate command can be omitted if the file is already name ordered or collated:
	# collate ensures that reads of the same name are grouped together in contiguous groups, 
	# but doesn't make any guarantees about the order of read names between groups
	samtools collate -o ${name}_collate.bam ${name}.bam
	# Add ms and MC tags for markdup to use later:
	# samtools fixmate – fills in mate coordinates and insert size fields
	samtools fixmate -m ${name}_collate.bam ${name}_fixmate.bam
	# Markdup needs position order:
	samtools sort -o ${name}_sort.bam ${name}_fixmate.bam
	# Finally mark duplicates:
	samtools markdup ${name}_sort.bam ${name}_markdup.bam

done


# Run samtools flagstat on the alignment file with marked duplicates
samtools flagstat ${dir_align}/son_markdup.bam

# to get the insert size average and std
samtools stats ${dir_align}/mother.sam | grep 'insert'

# For variant analysis, it’s important to know which read came from which sample.
# Therefore we add a tag to each read specifying its origin.

for person in 'father' 'mother' 'son' ; do
	samtools addreplacerg \
	-r ID:${person} \
	-r SM:${person} \
	-r PL:ILLUMINA \
	-o ${dir_align}/${person}_markdup.rg.bam \
	${dir_align}/${person}_markdup.bam
done


# view the header 
samtools view -H  son_markdup.rg.bam
# first few alignments
samtools view son_markdup.rg.bam | head


# An indexing can be compared to a kind of ‘phonebook’ of your sequence alignment file. 
# Indexing can be done with samtools as well, but it first needs to be sorted on coordinate
for person in 'father' 'mother' 'son' ; do
	samtools index ${dir_align}/${person}_markdup.bam
done











