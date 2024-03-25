#!/bin/bash
#SBATCH --job-name=HW_4		                                        # Job name
#SBATCH --partition=batch		                                    # Partition (queue) name
#SBATCH --ntasks=1			                                        # Single task job
#SBATCH --cpus-per-task=10		                                    # Number of cores per task 
#SBATCH --mem=40gb			                                        # Total memory for job
#SBATCH --time=2:00:00  		                                    # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/Homework_4/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/Homework_4/log.%j.err        # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                                 # Where to send mail 
#SBATCH --mail-type=END,FAIL                                        # Mail events (BEGIN, END, FAIL, ALL)

# set output directory variable
OUTDIR="/work/gene8940/fg69001/Homework_4"

# Using the SRA-toolkit to download and validate public NGS data
# Load modules
module load SRA-Toolkit/3.0.1-centos_linux64
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.15.1-GCC-11.3.0

# download paired-end illumina reads from NCBI for Ecoli strain C600
prefetch -O /work/gene8940/fg69001/Homework_4 SRR8082143 

# validate the downloaded data
vdb-validate /work/gene8940/fg69001/Homework_4 SRR80821

# Extract the downloaded data in SRA format to fastq format
fastq-dump --split-files --gzip /work/gene8940/fg69001/Homework_4/SRR8082143 -O /work/gene8940/fg69001/Homework_4/

# download the unmasked release 58 version of the E.coli MG1655 genome sequence files from Ensembl
SEQ="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
curl -s $SEQ | gunzip -c > $OUTDIR/ecoli_MG1655_ver58.fa

# construct reference genome index file
bwa index $OUTDIR/ecoli_MG1655_ver58.fa

# map paired-end sequencing reads to the reference genome
bwa mem -t 10 $OUTDIR/ecoli_MG1655_ver58.fa $OUTDIR/SRR8082143_1.fastq.gz $OUTDIR/SRR8082143_2.fastq.gz > $OUTDIR/SRR8082143.sam

# manipulate mapped reads using SAMtools
# convert SAM format to BAM format
samtools view $OUTDIR/SRR8082143.sam -O BAM -o $OUTDIR/SRR8082143.bam

# sort the BAM file
samtools sort --threads 6 $OUTDIR/SRR8082143.bam -o $OUTDIR/SRR8082143.sorted.bam

# Making BAM index from sorted BAM file of mapped reads
samtools index -@ 10 $OUTDIR/SRR8082143.sorted.bam

# Variant calling from short-read data

# call genotype likelihoods with `bcftools mpileup` using high mapping quality reads and create mpileup in gzipped VCF format
bcftools mpileup -Oz --threads 10 --min-MQ 60 -f $OUTDIR/ecoli_MG1655_ver58.fa $OUTDIR/SRR8082143.sorted.bam > $OUTDIR/SRR8082143.sorted.mpileup.vcf.gz

# call variants
bcftools call -Oz -m -v --threads 10 --ploidy 1 $OUTDIR/SRR8082143.sorted.mpileup.vcf.gz > $OUTDIR/SRR8082143.sorted.mpileup.call.vcf.gz

# filter variants with quality score less than 40 and those supported by fewer than 10 high quality mapped reads
bcftools filter -Oz -e 'QUAL<40 || DP<10' $OUTDIR/SRR8082143.sorted.mpileup.call.vcf.gz > $OUTDIR/SRR8082143.sorted.mpileup.call.filter.vcf.gz

# Generate an IGV readable index file
bcftools index $OUTDIR/SRR8082143.sorted.mpileup.call.filter.vcf.gz