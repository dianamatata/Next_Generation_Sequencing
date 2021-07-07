#!/bin/bash

### FOLDERS and FILES
home_dir=/home/davalos/workdir
dir_var=$home_dir/$dir_var

### Hard filering

# The developers of gatk strongly advise to do Variant Quality Score Recalibration (VQSR) for filtering SNPs and INDELs.
# VQSR stands for Variant Quality Score Recalibration. In a nutshell, it is a sophisticated filtering technique applied on the variant callset that uses machine learning to model the technical profile of variants in a training set and uses that to filter out probable artifacts from the callset.
# However, this is not always possible. For example, in the case of limited data availability and/or in the case you are working with non-model organisms and/or in the case you are a bit lazy and okay with a number of false positives.
# Our dataset is too small to apply VQSR. We will therefore do hard filtering instead.
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
# https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants

# Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION

for variant in 'SNP' 'INDEL' ; do
	gatk SelectVariants \
	--variant $dir_var/trio.vcf \
	--select-type-to-include $variant \
	--output $dir_var/trio.${variant}.vcf
done

## Filtering SNPs

# The command gatk VariantFiltration enables you to filter for both the INFO field (per variant) and FORMAT field (per genotype). For now we’re only interested in filtering variants. 
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

gatk VariantFiltration \
--variant $dir_var/trio.SNP.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output $dir_var/trio.SNP.filtered.vcf

## Filtering INDELs

gatk VariantFiltration \
--variant $dir_var/trio.INDEL.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--output $dir_var/trio.INDEL.filtered.vcf


## QualByDepth (QD)
# This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.
# Notice the values can be anywhere from 0-40

## FisherStrand (FS)
# This is the Phred-scaled probability that there is strand bias at the site. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. When there little to no strand bias at the site, the FS value will be close to 0.

## StrandOddsRatio (SOR)
# This is another way to estimate strand bias using a test similar to the symmetric odds ratio test. SOR was created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction and FS gives those variants a bad score. SOR will take into account the ratios of reads that cover both alleles.

## RMSMappingQuality (MQ)
# This is the root mean square mapping quality over all the reads at the site. Instead of the average mapping quality of the site, this annotation gives the square root of the average of the squares of the mapping qualities at the site. It is meant to include the standard deviation of the mapping qualities. Including the standard deviation allows us to include the variation in the dataset. A low standard deviation means the values are all close to the mean, whereas a high standard deviation means the values are all far from the mean.When the mapping qualities are good at a site, the MQ will be around 60.

## MappingQualityRankSumTest (MQRankSum)
# This is the u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping qualities of the reads supporting the reference allele and the alternate allele. A positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele; a negative value indicates the mapping qualities of the reference allele are higher than those supporting the alternate allele. A value close to zero is best and indicates little difference between the mapping qualities.

# no diff btw vcf and filtered.vcf in the number of records (lines) but we can check the filtered records (PASS) with:
grep -v "^#" variants/trio.SNP.filtered.vcf | cut -f 7 | sort | uniq -c
cat trio.INDEL.filtered.vcf | grep -v '^#' | cut -f7 | sort | uniq -c


## Merging filtered SNPs and INDELs

gatk MergeVcfs \
--INPUT $dir_var/trio.SNP.filtered.vcf \
--INPUT $dir_var/trio.INDEL.filtered.vcf \
--OUTPUT $dir_var/trio.filtered.vcf


### Evaluation by concordance

# For this region we have a highly curated truth set for the mother available. It originates from the Illumina Platinum truth set. You can find it at data/variants/NA12878.vcf.gz
# To check how well we did, we’d first need to extract a vcf with only the information of the mother.

# First select only SNPs and INDELs from the mother from the unfiltered vcf

gatk SelectVariants \
--variant $dir_var/trio.filtered.vcf \
--sample-name mother \
--exclude-non-variants \
--remove-unused-alternates \
--output $dir_var/mother.trio.filtered.vcf


cat $dir_var/mother.trio.filtered.vcf | grep -v '^#' | cut -f7 | sort | uniq -c

# Get the concordance with the truth set:

gatk Concordance \
--evaluation $dir_var/mother.trio.filtered.vcf \
--truth $home_dir/data/variants/NA12878.vcf.gz \
--intervals chr20:10018000-10220000 \
--summary $dir_var/concordance.mother.trio

cat $dir_var/concordance.mother.trio 
: <<'END'
type    TP      FP      FN      RECALL  PRECISION
SNP     319     5       9       0.973   0.985
INDEL   63      20      6       0.913   0.759
END

# 9 false negatives, i.e. SNPs we didn’t detect.
# https://en.wikipedia.org/wiki/Precision_and_recall
# thanks to filtering we removed 2FP. 


