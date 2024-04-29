library(DESeq2)
library(ggplot2)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(EnhancedVolcano)

#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results

countData_cotton <- read.csv("gossypium_bowtie_cotton_cvsn.matrix", header=T, row.names=1, sep="\t")

dim(countData_cotton)

head(countData_cotton)

barplot(colSums(countData_cotton)*1e-6, names=colnames(countData_cotton), ylab="Library size (millions)")

colData_cotton <- DataFrame(condition=factor(c("CCI","CCI","CCI","CCI","CNI","CNI","CNI")))

dds_cotton <- DESeqDataSetFromMatrix(countData=round(countData_cotton), colData_cotton, formula(~ condition))

#next two steps are for visualization of libraries

vsd_cotton <- vst(dds_cotton, blind = FALSE)
labels_cotton <- paste(c("CCI1", "CCI21", "CCI3", "CNI1", "CNI11", "CNI2"))
pc_cotton <- plotPCA(vsd_cotton, intgroup=c("condition")) + geom_text(label=labels_cotton)
gridExtra::grid.arrange(egg::set_panel_size(p=pc_cotton, width=unit(20, "cm"), height=unit(10, "cm")))
ggsave(filename = "pc12_cotton_cvsn.pdf", plot = egg::set_panel_size(p=pc_cotton, width=unit(12, "cm"), height=unit(10, "cm")))

#Next steps are for differential expression analysis in form of treated val untreated

#fold change log2(treated/untreated)

dds_cotton$condition <- relevel(dds_cotton$condition, ref = "CNI")

dds_cotton <- DESeq(dds_cotton)

par(mar=c(8,4,1,2))

boxplot(log10(assays(dds_cotton)[["cooks"]]), range=0, las=2, xlab = "", xaxt = "n")
axis(1, labels = FALSE)
text(x =  seq_along(labels_cotton), y = par("usr")[3] - 1, srt = 45, adj = 1,
     labels = labels_cotton, xpd = TRUE)

resultsNames(dds_cotton)

#this gives log2(n+1)
ntd_cotton <- normTransform(dds_cotton)

vsd_cotton <- vst(dds_cotton, blind=FALSE)
rld_cotton <- rlog(dds_cotton, blind=FALSE)
meanSdPlot(assay(rld_cotton))

#heatmap of count matrix
#select_val <- order(rowMeans(counts(dds_val,normalized=TRUE)),
#                    decreasing=TRUE)[1:20]
#df_val <- as.data.frame(colData(dds_val)[,c("condition")])
#pheatmap(assay(dds_val)[select_val,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df_val)

#heatmap of relatedness within samples
sampleDists_cotton <- dist(t(assay(vsd_cotton)))
sampleDistMatrix_cotton <- as.matrix(sampleDists_cotton)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

rownames(sampleDistMatrix_cotton) <- paste(c("CCI1","CCI21","CCI3","CCI4","CNI1","CNI11","CNI2"))
colnames(sampleDistMatrix_cotton) <- NULL
pheatmap(sampleDistMatrix_cotton,
         clustering_distance_rows=sampleDists_cotton,
         clustering_distance_cols=sampleDists_cotton,
         col=colors)

#plot fold counts across groups
res_cotton <- results(dds_cotton, contrast=c("condition","CCI","CNI"), alpha=0.001)

resLFC_cotton <- lfcShrink(dds_cotton, coef="condition_CCI_vs_CNI", type="apeglm")

plotCounts(dds_cotton, gene=which.min(res_cotton$padj), intgroup="condition")

#visualize deferentially expressed genes
MA_cotton <- plotMA(res_cotton)
MA_2_cotton <- plotMA(resLFC_cotton)
#gridExtra::grid.arrange(egg::set_panel_size(p=MA_val, width=unit(20, "cm"), height=unit(10, "cm")))
#ggsave(filename = "MA_cotton.png", plot = egg::set_panel_size(p=MA_val, width=unit(12, "cm"), height=unit(10, "cm")))
#order by adjusted p-value

resOrdered_cotton <- res_cotton[order(res_cotton$padj),]

head(resOrdered_cotton)

#get deferentially expressed genes matrix at FDR=5%, fold-change>2

sig_cotton <- resOrdered_cotton [!is.na(resOrdered_cotton$padj) & resOrdered_cotton$padj<0.05 & abs(resOrdered_cotton$log2FoldChange)>=0,]

head(sig_cotton)
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
                      legendIconSize =0,
                      col = c("grey30", "forestgreen", "royalblue", "red2"))

p1
ggsave(filename = "volcanoplot_cotton_cvsn.pdf", plot = egg::set_panel_size(p=p1, width=unit(10, "cm"), height=unit(8, "cm")))


write.csv(sig_cotton,file="cotton_cvsn_full.csv")

save.image("DESeq_cotton_3.23.22.RData")
