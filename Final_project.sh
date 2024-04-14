#!/bin/bash
#SBATCH --job-name=Final_project		                            # Job name
#SBATCH --partition=batch		                                    # Partition (queue) name
#SBATCH --ntasks=1			                                        # Single task job
#SBATCH --cpus-per-task=10		                                    # Number of cores per task 
#SBATCH --mem=40gb			                                        # Total memory for job
#SBATCH --time=2:00:00  		                                    # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/Final_project/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/Final_project/log.%j.err        # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                                 # Where to send mail 
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

# Load modules
module load SRA-Toolkit/3.0.1-centos_linux64
module load kallisto/0.48.0-gompi-2022a

OUTDIR="/work/gene8940/fg69001/Final_project/kallisto"
mkdir -p $OUTDIR

# generate a transcript De Bruijn Graph (T-DBG) stored in an index file
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_rna.fna.gz"
curl -s $URL | gunzip -c > $OUTDIR/cotton_rna.fa
kallisto index -i $OUTDIR/cotton_rna.fa.idx $OUTDIR/cotton_rna.fa

# Perform pseudoalignment and transcript abundance estimation using Kallisto T-DBG
THREADS=6
for i in CCI1, CCI21, CCI3, CNI1, CNI2, CNI11
do
  kallisto quant -t $THREADS -b 100 -i $OUTDIR/cotton_rna.fa.idx -o $OUTDIR/$i /work/gene8940/fg69001/Final_project/cotton_sequences/${i}_1.fq.gz /work/gene8940/fg69001/Final_project/cotton_sequences/${i}_2.fq.gz
done

# Create and activate a conda environment to perform differential expression using Sleuth
conda update -y conda
conda create -y --name sleuth_project
source activate sleuth_project
conda install -y --channel bioconda r-sleuth

source activate sleuth_project
R --no-save < /home/fg69001/GENE8940/Final_project.r