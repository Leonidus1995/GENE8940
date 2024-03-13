#!/bin/bash
#SBATCH --job-name=HW_3		                                        # Job name
#SBATCH --partition=batch		                                    # Partition (queue) name
#SBATCH --ntasks=1			                                        # Single task job
#SBATCH --cpus-per-task=10		                                    # Number of cores per task 
#SBATCH --mem=40gb			                                        # Total memory for job
#SBATCH --time=2:00:00  		                                    # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/fg69001/Homework_3/log.%j.out		# Location of standard output log files 
#SBATCH --error=/work/gene8940/fg69001/Homework_3/log.%j.err        # Location of standard error log files
#SBATCH --mail-user=fg69001@uga.edu                                 # Where to send mail 
#SBATCH --mail-type=END,FAIL                                        # Mail events (BEGIN, END, FAIL, ALL)

# set output directory variable
OUTDIR="/work/gene8940/fg69001/Homework_3" 

# Load the modules
module load canu/2.2-GCCcore-11.2.0
module load SPAdes/3.15.5-GCC-11.3.0
module load QUAST/5.2.0-foss-2022a
module load MUMmer/4.0.0rc1-GCCcore-11.3.0

# Long-read whole genome assembly with Canu
pacbio_file="/work/gene8940/instructor_data/ecoli_p6_25x.filtered.fastq.gz"
canu -p ecoli -d $OUTDIR/canu_ecoli_output_dir genomeSize=4.8m useGrid=false -pacbio-raw $pacbio_file

# Short-read whole genome assembly with SPAdes
illumina_file_1="/work/gene8940/instructor_data/s_6_1.fastq.gz"
illumina_file_2="/work/gene8940/instructor_data/s_6_2.fastq.gz"
spades.py -t 10 -k 21,33,55,77 --isolate --memory 40 --pe1-1 $illumina_file_1 --pe1-2 $illumina_file_2 -o $OUTDIR/spades_ecoli_output_dir

# Assembly quality assessment statistics using QUAST
SEQ="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-58/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz"
curl -s $SEQ | gunzip -c > $OUTDIR/ecoli_MG1655_ver58.fa
quast.py -o $OUTDIR/quast_output_dir -t 10 -r $OUTDIR/ecoli_MG1655_ver58.fa $OUTDIR/canu_ecoli_output_dir/ecoli.contigs.fasta $OUTDIR/spades_ecoli_output_dir/scaffolds.fasta

# Generate mummerplots for PacBio and Illumina assemblies
# For Illumina assembly
mkdir $OUTDIR/mummer
nucmer -t 10 $OUTDIR/ecoli_MG1655_ver58.fa $OUTDIR/spades_ecoli_output_dir/scaffolds.fasta -p $OUTDIR/mummer/mummer_ecoli_spades
delta-filter -1 $OUTDIR/mummer/mummer_ecoli_spades.delta > $OUTDIR/mummer/mummer_ecoli_spades.1delta
show-coords $OUTDIR/mummer/mummer_ecoli_spades.1delta
mummerplot --size large -layout --color -f --png $OUTDIR/mummer/mummer_ecoli_spades.1delta -p $OUTDIR/mummer/mummer_ecoli_spades

# For PacBio assembly
nucmer -t 10 $OUTDIR/ecoli_MG1655_ver58.fa $OUTDIR/canu_ecoli_output_dir/ecoli.contigs.fasta -p $OUTDIR/mummer/mummer_ecoli_canu
delta-filter -1 $OUTDIR/mummer/mummer_ecoli_canu.delta > $OUTDIR/mummer/mummer_ecoli_canu.1delta
show-coords $OUTDIR/mummer/mummer_ecoli_canu.1delta
mummerplot --size large -layout --color -f --png $OUTDIR/mummer/mummer_ecoli_canu.1delta -p $OUTDIR/mummer/mummer_ecoli_canu

