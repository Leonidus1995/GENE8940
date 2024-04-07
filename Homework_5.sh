#!/bin/bash
#SBATCH --job-name=HW_5		                                        # Job name
#SBATCH --partition=batch		                                    # Partition (queue) name
#SBATCH --ntasks=1			                                        # Single task job
#SBATCH --cpus-per-task=10		                                    # Number of cores per task 
#SBATCH --mem=40gb			                                        # Total memory for job
#SBATCH --time=2:00:00  		                                    # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/Homework_5/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/Homework_5/log.%j.err        # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                                 # Where to send mail 
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

# Load modules
module load SRA-Toolkit/3.0.1-centos_linux64
module load kallisto/0.48.0-gompi-2022a


OUTDIR="/work/gene8940/fg69001/Homework_5/kallisto"
mkdir -p $OUTDIR

# generate a transcript De Bruijn Graph (T-DBG) stored in an index file
curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz > $OUTDIR/ecoli_MG1655_cds.fa
kallisto index -i $OUTDIR/ecoli_MG1655_cds.fa.idx $OUTDIR/ecoli_MG1655_cds.fa

# Perform pseudoalignment and transcript abundance estimation using Kallisto T-DBG
THREADS=6
for i in SRR5344681 SRR5344682 SRR5344683 SRR5344684
do
  kallisto quant -t $THREADS -b 100 -i $OUTDIR/ecoli_MG1655_cds.fa.idx -o $i /work/gene8940/instructor_data/${i}_1.fastq.gz /work/gene8940/instructor_data/${i}_2.fastq.gz
done


# Create and activate a conda environment to perform differential expression using Sleuth
conda update -y conda
conda create -y --name sleuth_project
source activate sleuth_project
conda install -y --channel bioconda r-sleuth

source activate sleuth_project
R --no-save < /home/fg69001/GENE8940/homework5.r



