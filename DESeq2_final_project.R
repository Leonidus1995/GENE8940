library(rhdf5)
library(tximportData)
library(tximport)
library(readr)
library(biomaRt)
library(DESeq2)
library(PCAtools)
library(EnhancedVolcano)
library("pheatmap")

dir <- "/work/gene8940/fg69001/Final_project/kallisto"
setwd(dir)

samples <- read.csv(file.path(dir, "samples.csv"), header = TRUE)

files <- file.path(dir, "kallisto", samples$SampleName, "abundance.tsv")
names(files) <- samples$SampleName

# Mapping transcripts to genes
txdb <- makeTxDbFromGFF("")

