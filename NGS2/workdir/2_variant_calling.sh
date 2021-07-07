#!/bin/bash

### FOLDERS and FILES

home_dir=/home/davalos/workdir
dir_data=$home_dir/data
dir_var=$dir_data/variants
dir_align=$home_dir/alignment_pipe
ref_file=$dir_data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa

# vcfs:
        # A part of the dbsnp database: variants/GCF.38.filtered.renamed.vcf
        # A part of the 1000 genomes indel golden standard: variants/1000g_gold_standard.indels.filtered.vcf


### INDEXING

# Many algorithms work faster, or only work with an index of their input files
# samtools faidx – indexes or queries regions from a fasta file
# gatk CreateSequenceDictionary Creates a sequence dictionary for a reference sequence

samtools faidx $ref_file
gatk CreateSequenceDictionary --REFERENCE $ref_file

# input vcf files need to be indexed

for variant_file in $dir_var/*.vcf ; do
	gatk IndexFeatureFile --input $variant_file
done


### UNIFORM NAMING

# Before you start the alignment, it’s wise to check out what chromosome naming your input files are using in .fasta, .bam and .vcf files
# can be named as: 20, chr20, NC_000020.11 (Refseq name)
# changing chromosome names in a .fasta file is easier than in a .bam file
# in vcf: bcftools annotate --rename-chrs <tab-delimited-renaming> <input.vcf>
# in fasta: sed s/^>/>chr/g <reference.fasta>


### Base Quality Score Recalibration (BQSR)

# BQSR evaluates the base qualities on systematic error. It can ignore sites with known variants. BQSR helps to identify faulty base calls, and therefore reduces the chance on discovering false positive variant positions. Done in 2 steps

# 1 Recalibration with gatk BaseRecalibrator
# https://gatk.broadinstitute.org/hc/en-us/articles/360037593511-BaseRecalibrator
# First pass of the base quality score bqsr. Generates a bqsr table based on various covariates. The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

# 2 output application to the bam file with gatk ApplyBQSR
# https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR
# This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the base qualities of the input reads based on the bqsr table produced by the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

mkdir ${home_dir}/bqsr

for sample in 'father' 'mother' 'son' ; do

	gatk BaseRecalibrator \
	-I ${dir_align}/${sample}.bam \ 
	-R $ref_file \
	--known-sites ${dir_var}/1000g_gold_standard.indels.filtered.vcf \
	--known-sites ${dir_var}/GCF.38.filtered.renamed.vcf \
	-O ${home_dir}/bqsr/${sample}.recal.table

	gatk ApplyBQSR \
	-R $ref_file \
	-I ${dir_align}/${sample}.bam \
	--bqsr-recal-file ${home_dir}/bqsr/${sample}.recal.table \
	-O ${home_dir}/bqsr/${sample}.recal.bam

done


### Variant calling

# The command gatk HaplotypeCaller performs the actual variant calling in gatk
# HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
# https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#  In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper
# We’ll focus on a small region, so add --intervals chr20:10018000-10220000

dir_out_var=/home/davalos/workdir/variants
mkdir $dir_out_var

for sample in 'father' 'mother' 'son' ; do
		
	gatk HaplotypeCaller \
	   -R $ref_file \
	   -I ${home_dir}/bqsr/${sample}.recal.bam \
	   --intervals chr20:10018000-10220000 \
	   -O $dir_out_var/${sample}.HC.vcf.gz

done

#  get the number of records in a vcf
for sample in 'father' 'mother' 'son' ; do
	echo $sample
	gunzip -d $dir_out_var/${sample}.HC.vcf.gz
	grep -v '^#' $dir_out_var/${sample}.HC.vcf | wc -l
	# bgzip $dir_out_var/${sample}.HC.vcf
done
# father463 mother411 son421


# You can get some more statistics with gatk VariantsToTable

for sample in 'father' 'mother' 'son' ; do
	gatk VariantsToTable \
	--variant $dir_out_var/${sample}.HC.vcf \
	--fields CHROM -F POS -F TYPE -GF GT \
	--output $dir_out_var/${sample}.HC.table
done


# get the number of SNPs and INDELs
grep -c "SNP" $dir_out_var/${sample}.HC.table
grep -c "INDEL" $dir_out_var/${sample}.HC.table
cut -f 3 $dir_out_var/${sample}.HC.table | tail -n +2 | sort | uniq -c


# call again gatk HaplotypeCaller but generating GVCFs
# In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).

dir_GVCF=/home/davalos/workdir/variants
mkdir $dir_GVCF

for sample in 'father' 'mother' 'son' ; do
        gatk HaplotypeCaller  \
        --reference $ref_file \
        --input ${home_dir}/bqsr/${sample}.recal.bam \
        --output $dir_out_var/${sample}.HC.g.vcf \
        --bam-output $dir_out_var/${sample}.phased.bam \
        --intervals chr20:10018000-10220000 \
        --emit-ref-confidence GVCF
done


### Combining GVCFs

# Now that we have all three GVCFs of the mother, father and son, we can combine them into a database. 
# We will do the variant calling on all three samples. Later we want to combine the variant calls. For efficient merging of vcfs, we will need to output the variants as a GVCF. To do that, we will use the option --emit-ref-confidence GVCF. Also, we’ll visualise the haplotype phasing with IGV in the next section. For that we’ll need a phased bam. You can get this output with the argument --bam-output.

gatk GenomicsDBImport \
--variant  $dir_out_var/mother.HC.g.vcf \
--variant  $dir_out_var/father.HC.g.vcf \
--variant  $dir_out_var/son.HC.g.vcf \
--intervals chr20:10018000-10220000 \
--genomicsdb-workspace-path genomicsdb


# You can retrieve the combined vcf from the database with gatk GenotypeGVCFs

gatk GenotypeGVCFs \
--reference $ref_file \
--variant gendb://$home_dir/genomicsdb \
--intervals chr20:10018000-10220000 \
--output $dir_out_var/trio.vcf


