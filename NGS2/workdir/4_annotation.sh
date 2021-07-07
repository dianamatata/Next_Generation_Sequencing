#!/bin/bash

### FOLDERS and FILES
# To use the human genome is a reference, we have downloaded the database with:
# snpEff download -v GRCh38.99

mkdir annotation

snpEff -Xmx4g \
	-v \
	-o gatk \
	GRCh38.99 \
	variants/trio.filtered.vcf > annotation/trio.filtered.snpeff.vcf


