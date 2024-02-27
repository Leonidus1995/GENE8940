#!/bin/bash
#SBATCH --job-name=HW_2		                            # Job name
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
OUTDIR="/work/gene8940/fg69001/Homework_2" 
