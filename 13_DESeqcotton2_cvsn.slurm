#!/bin/bash
#SBATCH --job-name=7_DESeq2_val_to_val
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --time=48:00:00

cd $SLURM_SUBMIT_DIR
module load R/4.1.0-foss-2019b

R CMD BATCH 12_DESeqcotton_cvsn.R
