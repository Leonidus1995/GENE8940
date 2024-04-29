#!/bin/bash
#SBATCH --job-name=Bowtie		                            
#SBATCH --partition=batch		                                    
#SBATCH --ntasks=1			                                        
#SBATCH --cpus-per-task=10		                                     
#SBATCH --mem=40gb			                                        
#SBATCH --time=2:00:00  		                                    
#SBATCH --output=/work/gene8940/fg69001/Final_project/log.%j.out		 
#SBATCH --error=/work/gene8940/fg69001/Final_project/log.%j.err        
#SBATCH --mail-user=fg69001@uga.edu                                  
#SBATCH --mail-type=END,FAIL  


module load SRA-Toolkit/3.0.1-centos_linux64
module load Bowtie2/2.5.2-GCC-11.3.0
module load SAMtools/1.14-GCC-11.2.0

bowtie2-build /work/gene8940/fg69001/Final_project/bowtie/cotton_rna.fa /work/gene8940/fg69001/Final_project/bowtie/ref_index/cotton 

reference_index="/work/gene8940/fg69001/final_project_test/reference_index"
for i in CCI1 CCI21 CCI3 CNI1 CNI2 CNI11
do
    bowtie2 -x $reference_index/cotton_index.idx -1 /work/gene8940/fg69001/Final_project/cotton_sequences/${i}_1.fq.gz -2 /work/gene8940/fg69001/Final_project/cotton_sequences/${i}_2.fq.gz -S /work/gene8940/fg69001/Final_project/bowtie/mapped/${i}.sam --threads 10
done


for i in CCI1 CCI21 CCI3 CNI1 CNI2 CNI11
do
    samtools view -b -o /work/gene8940/fg69001/Final_project/bowtie/mapped/${i}.bam /work/gene8940/fg69001/Final_project/bowtie/mapped/${i}.sam
done