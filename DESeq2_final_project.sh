#!/bin/bash
#SBATCH --job-name=DESeq2_final_project		                            
#SBATCH --partition=batch		                                    
#SBATCH --ntasks=1			                                        
#SBATCH --cpus-per-task=24		                                     
#SBATCH --mem=50gb			                                        
#SBATCH --time=72:00:00  		                                    
#SBATCH --output=/work/gene8940/fg69001/Final_project/log.%j.out		 
#SBATCH --error=/work/gene8940/fg69001/Final_project/log.%j.err        
#SBATCH --mail-user=fg69001@uga.edu                                  
#SBATCH --mail-type=END,FAIL 

module load R/4.3.2-foss-2022b

R CMD BATCH /home/fg69001/GENE8940/DESeq2_final_project.R
