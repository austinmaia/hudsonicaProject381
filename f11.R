setwd("~/Documents/GitHub/hudsonicaProject381/data")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(hexbin)
library(eulerr)
library(pheatmap)

countsTable <- read.table("salmon.gene.counts.f11.txt", header=TRUE, row.names=1)
head(countsTable)

dim(countsTable) # 130580     6

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

conds <- read.delim("F11.txt", header=T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ treatment)
dds <- dds[rowSums(counts(dds)) > 60]
dds <- DESeq(dds)
vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("treatment"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Differential expression
resultsNames(dds)

#Combined vs ambient
resF11_Combo <- results(dds, alpha=0.05)
resF11_Combo <- resF11_Combo[order(resF11_Combo$padj),] #ordering rows by adjusted p-value
head(resF11_Combo)
# log2 fold change (MLE): treatment OWA vs AM 
# Wald test p-value: treatment OWA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN42020_c1_g1   521.568        3.89502  0.460805   8.45265
# TRINITY_DN39749_c0_g1   175.626      -11.07517  1.390028  -7.96759
# TRINITY_DN30572_c0_g1   155.182       10.63355  1.357951   7.83058
# TRINITY_DN13962_c0_g2   163.068       10.70457  1.378184   7.76716
# TRINITY_DN35884_c0_g1   162.992       -2.45551  0.320239  -7.66776
# TRINITY_DN34691_c0_g1  2271.224        1.80119  0.235867   7.63645
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN42020_c1_g1 2.84772e-17 2.12745e-12
# TRINITY_DN39749_c0_g1 1.61799e-15 6.04375e-11
# TRINITY_DN30572_c0_g1 4.85611e-15 1.20928e-10
# TRINITY_DN13962_c0_g2 8.02683e-15 1.49915e-10
# TRINITY_DN35884_c0_g1 1.75026e-14 2.61513e-10
# TRINITY_DN34691_c0_g1 2.23297e-14 2.78031e-10

summary(resF11_Combo)
# out of 75787 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 839, 1.1%
# LFC < 0 (down)     : 1300, 1.7%
# outliers [1]       : 1080, 1.4%
# low counts [2]     : 0, 0%
# (mean count < 8)


resF11_Combo <- resF11_Combo[!is.na(resF11_Combo$padj),]
degsCombo <- row.names(resF11_Combo[resF11_Combo$padj < 0.05,])
length(degsCombo) 
#2139


#Heatmaps

topgenes <- head(rownames(resF11_Combo),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

fit1 <- euler(c("Combined" = 2139))

plot(fit1,  lty = 1:3, quantities = TRUE)
