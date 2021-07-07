# https://github.com/BrooksLabUCSC/flair

# FLAIR (Full-Length Alternative Isoform analysis of RNA) for the correction, isoform definition, and alternative splicing analysis of noisy reads. FLAIR has primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads.

alignment_folder=~/longreads2/aligned_results
SAMPLE1=ERR3577102
SAMPLE2=ERR3577098

# python flair.py align -g genome.fa -r <reads.fq>|<reads.fa> [options]
# python flair.py correct -q query.bed12 -g genome.fa [options] with gtf file in options
# python flair.py collapse -g genome.fa -r <reads.fq>|<reads.fa> -q <query.psl>|<query.bed> [options]
# python flair.py quantify -r reads_manifest.tsv -i isoforms.fasta [options]


for alignment_file in $alignment_folder/ERR3577098_G500k $alignment_folder/ERR3577102_G500k ; do
	python3 ~/flair/bin/bam2Bed12.py -i $alignment_file.bam > $alignment_file.bed12
	python3 ~/flair/bin/bam2Bed12.py align -g ~/longreads2/refs/Homo_sapiens.GRCh38.dna_rm.chromosome.12.fa.gz -r ~/longreads2/ERR3577098.fastq > ERR3577098_flair_aligned.sam

done
