#!/bin/bash
#SBATCH --job-name=1_rsem_3
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --mail-user=sp36142@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR
module load Bowtie2/2.4.1-GCC-8.3.0
module load RSEM/1.3.3-foss-2019b

#/apps/eb/RSEM/1.3.3-foss-2019b/bin/rsem-generate-data-matrix \
#gossypium_bowtie_CCI1.genes.results gossypium_bowtie_CCI21.genes.results gossypium_bowtie_CCI3.genes.results gossypium_bowtie_CCI4.genes.results \
#gossypium_bowtie_CMI11.genes.results gossypium_bowtie_CMI21.genes.results gossypium_bowtie_CMI31.genes.results gossypium_bowtie_CMI4.genes.results  \
#> gossypium_bowtie_cotton_cvsm.matrix

/apps/eb/RSEM/1.3.3-foss-2019b/bin/rsem-generate-data-matrix \
gossypium_bowtie_CCI1.genes.results gossypium_bowtie_CCI21.genes.results gossypium_bowtie_CCI3.genes.results gossypium_bowtie_CCI4.genes.results \
gossypium_bowtie_CNI1.genes.results gossypium_bowtie_CNI11.genes.results gossypium_bowtie_CNI2.genes.results   \
> gossypium_bowtie_cotton_cvsn.matrix

#/apps/eb/RSEM/1.3.3-foss-2019b/bin/rsem-generate-data-matrix \
#gossypium_bowtie_CMI11.genes.results gossypium_bowtie_CMI21.genes.results gossypium_bowtie_CMI31.genes.results gossypium_bowtie_CMI4.genes.results  \
#gossypium_bowtie_CNI1.genes.results gossypium_bowtie_CNI11.genes.results gossypium_bowtie_CNI2.genes.results   \
#> gossypium_bowtie_cotton_mvsn.matrix
