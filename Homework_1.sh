#!/bin/bash
#SBATCH --job-name=HW_1		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task 
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/log.%j.err       # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                     # Where to send mail 
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/fg69001/Homework_1"                  
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.58.gff3.gz"
curl -s $URL | gunzip -c > $OUTDIR/ecoli_MG1655.gff           

#if output directory doesn't exist, create it for sure
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

cut -f3 /work/gene8940/fg69001/Homework_1/ecoli_MG1655.gff | grep "CDS" | sort | uniq -c > /work/gene8940/fg69001/Homework_1/results.txt