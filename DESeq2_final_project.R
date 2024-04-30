library(tximport)
library(jsonlite)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(vsn)

dir <- "/work/gene8940/fg69001/Final_project/kallisto"
setwd(dir)

samples <- read.csv(file.path(dir, "samples.csv"), header = TRUE)

files <- file.path(dir, "kallisto", samples$SampleName, "abundance.tsv")
names(files) <- samples$SampleName

# Mapping transcripts to genes
txdb <- makeTxDbFromGFF(file.path(dir, "cotton.gtf.gz"))
saveDb(txdb, file.path(dir, "cotton.txdb"))

# Create tx2gene data frame
txdb <- loadDb(file.path(dir, "cotton.txdb"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Import Kallisto transcript quantifications into R,
# collapsing to the gene level using the information in tx2gene
txi <- tximport(files, type = "kallisto",
                tx2gene = tx2gene, ignoreAfterBar = TRUE)

# DESeq2 analysis
colData_cotton <- DataFrame(condition=factor(c("CCI","CCI","CCI","CNI","CNI","CNI")))
dds <- DESeqDataSetFromTximport(txi,
                                colData = colData_cotton,
                                formula(~ condition))


vsd <- vst(dds, blind = FALSE)

# Principal Component Plot
#pca_data <- plotPCA(vsd, intgroup = c("dex", "SampleName"), returnData = TRUE)
#percent_var <- round(100 * attr(pca_data, "percent_var"))
#ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = SampleName)) +
  #geom_point(size = 3) +
  #xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  #ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  #coord_fixed()

labels_cotton <- paste(c("CCI1", "CCI21", "CCI3", "CNI1", "CNI11", "CNI2"))
pc_cotton <- plotPCA(vsd, intgroup=c("condition")) + geom_text(label=labels_cotton)
gridExtra::grid.arrange(egg::set_panel_size(p=pc_cotton, width=unit(20, "cm"), height=unit(10, "cm")))
ggsave(filename = "pc12_cotton.pdf", plot = egg::set_panel_size(p=pc_cotton, width=unit(12, "cm"), height=unit(10, "cm")))

# Next steps are for differential expression analysis in form of treated val untreated

# fold change log2(treated/untreated)

dds$condition <- relevel(dds$condition, ref = "CNI")

dds_cotton <- DESeq(dds)

#this gives log2(n+1)
ntd_cotton <- normTransform(dds_cotton)

vsd_cotton <- vst(dds_cotton, blind = FALSE)



# plot fold counts across groups
res_cotton <- results(dds_cotton, contrast=c("condition","CCI","CNI"), alpha=0.001)

resLFC_cotton <- lfcShrink(dds_cotton, coef="condition_CCI_vs_CNI", type="apeglm")

plotCounts(dds_cotton, gene=which.min(res_cotton$padj), intgroup = "condition")

resOrdered_cotton <- res_cotton[order(res_cotton$padj), ]


#get deferentially expressed genes matrix at FDR=5%, fold-change>2

sig_cotton <- resOrdered_cotton[!is.na(resOrdered_cotton$padj) & resOrdered_cotton$padj<0.05 & abs(resOrdered_cotton$log2FoldChange)>=0, ]

#visualize volcano plot for differentially expressed genes

p1 <- EnhancedVolcano(res_cotton,
                      lab = rownames(res_cotton),
                      x = "log2FoldChange",
                      y = "pvalue",
                      pCutoff = 10e-4,
                      FCcutoff = 4,
                      title = c(),
                      subtitle = "Differential expression",
                      caption = c(), legendLabels = c(),
                      labSize = 0.0,
                      legendPosition = c(),
                      legendLabSize = 0,
                      legendIconSize = 0,
                      col = c("grey30", "forestgreen", "royalblue", "red2"))

p1
ggsave(filename = "volcanoplot_cotton.pdf", plot = egg::set_panel_size(p = p1, width = unit(10, "cm"), height = unit(8, "cm")))

write.csv(sig_cotton, file = "cotton_cvsn_full.csv")

save.image("DESeq_cotton_3.23.22.RData")