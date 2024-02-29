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

# Load modules to manipulate sequence and annotation files
module load BEDOPS/2.4.41-foss-2021b 
module load BEDTools/2.30.0-GCC-11.2.0
module load SAMtools/1.14-GCC-11.2.0
module load ucsc/443

# Convert the GFF file to BED format
convert2bed --input=gff < ecoli_MG1655_ver58.gff > $OUTDIR/ecoli58.bed

# Filter the BED file to create new BED file with only CDS regions
grep "ID=CDS" ecoli58.bed > $OUTDIR/ecoli58_1_cds.bed

# Create a genome index file for the genome sequence fasta file
samtools faidx ecoli_MG1655_ver58.fa

# Extract the first two columns of the index file and save it in the text file
cut -f1,2 ecoli_MG1655_ver58.fa.fai > $OUTDIR/ecoli_MG1655_index.txt

# Find the non CDS regions
bedtools complement -i ecoli58_1_cds.bed -g ecoli_MG1655_index.txt > $OUTDIR/ecoli_intergenic.bed

# Extract FASTA sequences for CDS and non-CDS regions from a reference genome based on respective annotated BED files
bedtools getfasta -fi ecoli_MG1655_ver58.fa -bed ecoli58_1_cds.bed -fo $OUTDIR/ecoli58_cds.fa
bedtools getfasta -fi ecoli_MG1655_ver58.fa -bed ecoli_intergenic.bed -fo $OUTDIR/ecoli58_intergenic.fa

# Calculating the %GC content of CDS and non-CDS regions
faCount ecoli58_cds.fa > $OUTDIR/results_cds.txt
awk 'NR > 1 { gc_count += $4 + $5; total_count += $6 + $3 + $4 + $5} END { 
    gc_content = (gc_count / total_count) * 100; 
    print("Overall GC content:", gc_content)}' results_cds.txt > $OUTDIR/summary_cds.txt

faCount ecoli58_intergenic.fa > $OUTDIR/results_intergenic.txt
awk 'NR > 1 { gc_count += $4 + $5; total_count += $6 + $3 + $4 + $5} END { 
    gc_content = (gc_count / total_count) * 100; 
    print("Overall GC content:", gc_content)}' results_intergenic.txt > $OUTDIR/summary_intergenic.txt

# faCount -summary command also give the similar results for GC content
faCount -summary ecoli58_cds.fa > $OUTDIR/facount_cds.txt
faCount -summary ecoli58_intergenic.fa > $OUTDIR/facount_intergenic.txt