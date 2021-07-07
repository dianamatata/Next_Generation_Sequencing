#!/bin/bash

ref=~/ecoli/ref_genome/ecoli-strK12-MG1655.fasta
output_name=~/ecoli/ref_genome/ecoli-strK12-MG1655.fasta
cmd="bowtie2-build $ref $output_name"
echo $cmd
eval $cmd
