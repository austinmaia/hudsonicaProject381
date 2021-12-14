
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

countsTable <- read.table("salmon.gene.counts.f4.txt", header=TRUE, row.names=1)
head(countsTable)

dim(countsTable) # 130580     12


countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

conds <- read.delim("ahud_samples_F4.txt", header=T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ treatment)
dds <- dds[rowSums(counts(dds)) > 120]
dds <- DESeq(dds)
vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("treatment"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()


#Combined vs ambient
resF4_Combo <- results(dds, alpha=0.05)
resF4_Combo <- resF4_Combo[order(resF4_Combo$padj),] #ordering rows by adjusted p-value
head(resF4_Combo)
# log2 fold change (MLE): treatment OWA vs AM 
# Wald test p-value: treatment OWA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN29519_c0_g1   55.5800       21.69690  2.325614   9.32953
# TRINITY_DN13334_c0_g1   93.4191       -3.85708  0.447504  -8.61910
# TRINITY_DN32344_c0_g1   27.3046       19.91202  2.384014   8.35231
# TRINITY_DN50669_c0_g1  110.4042       22.92111  2.849802   8.04305
# TRINITY_DN6992_c2_g1   118.5751       10.68015  1.367250   7.81141
# TRINITY_DN40540_c1_g1   35.3964      -22.45117  2.927096  -7.67012
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN29519_c0_g1 1.06338e-20 7.88509e-16
# TRINITY_DN13334_c0_g1 6.74838e-18 2.50200e-13
# TRINITY_DN32344_c0_g1 6.69427e-17 1.65462e-12
# TRINITY_DN50669_c0_g1 8.76273e-16 1.62441e-11
# TRINITY_DN6992_c2_g1  5.65525e-15 8.38685e-11
# TRINITY_DN40540_c1_g1 1.71838e-14 2.12366e-10


summary(resF4_Combo)
# out of 75477 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 305, 0.4%
# LFC < 0 (down)     : 335, 0.44%
# outliers [1]       : 1326, 1.8%
# low counts [2]     : 0, 0%
# (mean count < 9)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


#OA vs ambient
resF4_OA <- results(dds,contrast = c("treatment", "AM", "OA"), alpha=0.05)
resF4_OA <- resF4_OA[order(resF4_OA$padj),] #ordering rows by adjusted p-value
head(resF4_OA)
# log2 fold change (MLE): treatment AM vs OA 
# Wald test p-value: treatment AM vs OA 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN67785_c0_g1   39.5067      22.217584 2.4197646   9.18171
# TRINITY_DN46090_c4_g1   30.8544      21.299563 2.3717806   8.98041
# TRINITY_DN402_c0_g1   3519.7294      -0.428397 0.0488058  -8.77758
# TRINITY_DN29519_c0_g1   55.5800     -19.831984 2.3290044  -8.51522
# TRINITY_DN32344_c0_g1   27.3046     -19.722245 2.3850519  -8.26911
# TRINITY_DN11863_c1_g1   37.4748      22.639924 2.9904129   7.57084
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN67785_c0_g1 4.24288e-20 3.14614e-15
# TRINITY_DN46090_c4_g1 2.69758e-19 1.00014e-14
# TRINITY_DN402_c0_g1   1.67027e-18 4.12841e-14
# TRINITY_DN29519_c0_g1 1.66275e-17 3.08236e-13
# TRINITY_DN32344_c0_g1 1.34968e-16 2.00161e-12
# TRINITY_DN11863_c1_g1 3.70832e-14 4.58293e-10
summary(resF4_OA)
#out of 75477 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 197, 0.26%
# LFC < 0 (down)     : 316, 0.42%
# outliers [1]       : 1326, 1.8%
# low counts [2]     : 0, 0%
# (mean count < 9)


#OW vs ambient
resF4_OW <- results(dds,contrast = c("treatment", "AM", "OW"), alpha=0.05)
resF4_OW <- resF4_OW[order(resF4_OW$padj),] #ordering rows by adjusted p-value
head(resF4_OW)
# log2 fold change (MLE): treatment AM vs OW 
# Wald test p-value: treatment AM vs OW 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN34045_c0_g1   78.0749        23.6965  2.475966   9.57060
# TRINITY_DN29519_c0_g1   55.5800       -20.7972  2.326901  -8.93773
# TRINITY_DN32344_c0_g1   27.3046       -20.5371  2.382892  -8.61854
# TRINITY_DN480_c2_g2    150.4565        -2.7078  0.336422  -8.04881
# TRINITY_DN27511_c1_g1   58.5785        22.5659  2.920146   7.72766
# TRINITY_DN46926_c0_g1   40.2302        22.6365  2.937621   7.70571
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN34045_c0_g1 1.06282e-21 7.88095e-17
# TRINITY_DN29519_c0_g1 3.97251e-19 1.47283e-14
# TRINITY_DN32344_c0_g1 6.78136e-18 1.67615e-13
# TRINITY_DN480_c2_g2   8.36047e-16 1.54984e-11
# TRINITY_DN27511_c1_g1 1.09541e-14 1.60805e-10
# TRINITY_DN46926_c0_g1 1.30117e-14 1.60805e-10

summary(resF4_OW)
# out of 75477 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 182, 0.24%
# LFC < 0 (down)     : 259, 0.34%
# outliers [1]       : 1326, 1.8%
# low counts [2]     : 0, 0%



#Heatmaps

topgenes <- head(rownames(resF4_OA),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF4_Combo),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF4_OW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)


# Fixing venn diagram #
# rename padj columns
# merge data frames by trinity ID 
# filter NAs
# find degs

names(resF4_Combo)[6] <- 'pOWA'
names(resF4_OA)[6] <- 'pOA'
names(resF4_OW)[6] <- 'pOW'

fixed <- merge.data.frame(resF4_Combo, resF4_OA, by=0, all=TRUE)
fixed$pOW <- resF4_OW$pOW
row.names(fixed) <- fixed$Row.names
fixed <- fixed[!is.na(fixed$pOA),]
fixed <- fixed[!is.na(fixed$pOW),]
fixed <- fixed[!is.na(fixed$pOWA),]
degsCombo <- row.names(fixed[fixed$pOWA < 0.05,])
length(degsCombo) #638

degsOA <- row.names(fixed[fixed$pOA < 0.05,])
length(degsOA) #512

degsOW <- row.names(fixed[fixed$pOW < 0.05,])
length(degsOW) #431

## VENN DIAGRAM ##
# Total
length(degsOA)  # 512
length(degsOW)  # 431
length(degsCombo)  # 638

# Intersections
length(intersect(degsOA,degsOW))  # 5
length(intersect(degsOA,degsCombo))  # 160
length(intersect(degsCombo,degsOW))  # 5

intEL <- intersect(degsOA,degsOW)
length(intersect(degsCombo,intEL)) # 3

# Number unique
512-5-160-3 #OA 344
431-5-5-3 #OW 418
638-160-5-3 #Combo 470

# unique, present in 2, present in 3
fit1 <- euler(c("OA" = 344, "OW" = 418, "Combined" = 470, "OA&OW" = 5, "OA&Combined" = 160, "OW&Combined" = 5, "OA&OW&Combined" = 3))

plot(fit1,  lty = 1:3, quantities = TRUE)
