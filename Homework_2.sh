#!/bin/bash
#SBATCH --job-name=HW_2		                            # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task 
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/Homework_2/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/Homework_2/log.%j.err       # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                     # Where to send mail 
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# set output directory variable
OUTDIR="/work/gene8940/fg69001/Homework_2" 

SEQ="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
curl -s $SEQ | gunzip -c > $OUTDIR/ecoli_MG1655_ver58.fa
ANN="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.58.chromosome.Chromosome.gff3.gz"
curl -s $ANN | gunzip -c > $OUTDIR/ecoli_MG1655_ver58.gff

# Load modules to manipulate annotation file
module load BEDOPS/2.4.41-foss-2021b 
module load BEDTools/2.30.0-GCC-11.2.0
module load SAMtools/1.14-GCC-11.2.0
module load ucsc/443

# Convert the GFF file to BED format
convert2bed --input=gff < ecoli_MG1655_ver58.gff > ecoli58.bed

# Filter the BED file to create new BED file with only CDS regions
cut -f1 /work/gene8940/fg69001/Homework_2/ecoli58.bed | grep "CDS" | sort | uniq -c > /work/gene8940/fg69001/Homework_2/ecoli58_CDS.bed