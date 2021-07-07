ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.12.fa.gz### Aim: Align long reads from RNA-seq data to a reference genome.:

# In this project, you will be working with data from:
# Clark, M. B. et al (2020). Long-read sequencing reveals the complex splicing profile of the psychiatric risk gene CACNA1C in human brain. Molecular Psychiatry, 25(1), 37–47. https://doi.org/10.1038/s41380-019-0583-1.

### download data

# Download the human reference genome:

mkdir refs
cd refs
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
cd


#download samples from https://www.ncbi.nlm.nih.gov/bioproject/PRJEB34660

SAMPLE1=ERX3574810
SAMPLE2=ERX3574806

for sample in $SAMPLE1 $SAMPLE2 ; do
	prefetch $sample
	fastq-dump $sample
done
# note from Geert: remove --split-files from fastq-dump because not paired end reads because not relevant, we have single-end reads. but no need to rerun, the software can handle it
# note2: # download file: prefetch will download and save SRA file related to SRR accession in the current directory under newly created SRA accession directory NOTE: With fastq-dump and fasterq-dump, prefetch step is unncessary and you can directly download sequence data in FASTQ format


# perform fastqc on a *.fastq file with 2 different methods and output in 2 different folders
mkdir fastqc Nanoplot Nanoplot/ERR3577098 Nanoplot/ERR3577102
fastqc *.fastq -o fastqc/

SAMPLE1=ERR3577102
SAMPLE2=ERR3577098
for sample in $SAMPLE1 $SAMPLE2 ; do
	NanoPlot --fastq $sample.fastq -o Nanoplot/${sample}
done

### Analyze QC results
# Find the equation to calculate error probability from quality score
# https://en.wikipedia.org/wiki/Phred_quality_score


### Align with minimap2 with default parameters
minimap2 -ax map-ont Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ERX3574810_1.fastq > aligned_results/ERX3574810.sam
minimap2 -ax map-ont Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ERX3574806_1.fastq > aligned_results/ERX3574806.sam


### figure out parameters for minimap2
# Figure how you should set parameters -x and -G
# map-ont: Slightly more sensitive for Oxford Nanopore to reference mapping (-k15). For PacBio reads, HPC minimizers consistently leads to faster performance and more sensitive results in comparison to normal minimizers. For Oxford Nanopore data, normal minimizers are better, though not much. The effectiveness of HPC is determined by the sequencing error mode.
# -x for multiple parameters
# -a Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
# G: Maximum gap on the reference (effective with -xsplice/--splice). This option also changes the chaining and alignment band width to NUM. Increasing this option slows down spliced alignment. [200k]
# need to find the maximum intron size for CACNA1C in order to determine G.
# https://www.ncbi.nlm.nih.gov/gene?term=NM_001129830.3 page of the gene CACNA1C, longest intron 477 774 nucleotides

REFERENCE=refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
G_VALUE=500k
PARAMETER=map-ont
SAMPLE1=ERR3577102
SAMPLE2=ERR3577098

for sample in $SAMPLE1 $SAMPLE2 ; do
        fastqfile=${sample}.fastq
        echo $fastqfile

        minimap2 \
        -a \
        -x $PARAMETER \
        -G $G_VALUE \
        $REFERENCE \
        $fastqfile \
        | samtools sort \
        | samtools view -bh > ${sample}_G${G_VALUE}.bam

        samtools index ${sample}_G${G_VALUE}.bam
done

# samtools view *.bam, zcat would give a weird result
# get sam files
minimap2 -ax map-ont -G 500k $REFERENCE ERR3577098.fastq  | samtools sort > ERR3577098.sam
minimap2 -ax map-ont -G 500k $REFERENCE ERR3577102.fastq  | samtools sort > ERR3577102.sam
minimap2 -ax map-ont -G 500k $REFERENCE ERR3577098.fastq  | samtools sort | samtools view -bh > ERR3577098_500k_n.bam
minimap2 -ax map-ont -G 500k $REFERENCE ERR3577102.fastq  | samtools sort  | samtools view -bh > ERR3577102_500k_n.bam
samtools index ERR3577102_500k_n.bam
samtools index ERR3577098_500k_n.bam


# 500k understood, 500000 gives weird results.


# we have weird results, Geert version
#!/usr/bin/env bash
minimap2 \
-a \
-x splice \
-G 500k \
/home/lborcard/project_long_read/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
/home/lborcard/project_long_read/raw_data/ERR3577102.fastq \
| samtools sort - \
| samtools view -bh - \
> ~/CACNA1C/test.bam
samtools index ~/CACNA1C/test.bam

### end of Geert version

# analysis
samtools flagstat ERR3577098.sam > ERR3577098sam_flagstats.txt
samtools flagstat ERR3577102.sam > ERR3577102sam_flagstats.txt

samtools stats ERR3577098.sam  | grep ^SN | cut -f 2,3 > ERR3577098sam_stats.txt
samtools stats ERR3577102.sam  | grep ^SN | cut -f 2,3 > ERR3577102sam_stats.txt



# Reference genome (chromosome 12 only)
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.12.fa.gz
# GTF: 
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
# flair
# FLAIR (Full-Length Alternative Isoform analysis of RNA) for the correction, isoform definition, and alternative splicing analysis of noisy reads. FLAIR has primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads.
$alignment
python3 ~/flair/bin/bam2Bed12.py \
-i $alignment.bam \
> $alignment.bed12

