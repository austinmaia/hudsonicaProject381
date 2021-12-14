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

countsTable <- read.table("salmon.gene.counts.f04.txt", header=TRUE, row.names=1)
head(countsTable)

dim(countsTable) # 130580     24

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

conds <- read.delim("ahud_samples_f04.txt", header=T, stringsAsFactors = T, row.names = 1)
head(conds)
dim(conds)

dds <- DESeqDataSetFromMatrix(countData=countsTableRound, 
                              colData=conds, 
                              design=~treatment+generation
                              +treatment:generation)

### Filter out genes with too few reads - keep reads with average > 10 reads per sample

dds <- dds[rowSums(counts(dds)) > 240] #10 * 24 samples
dim(dds) #79036    24

dds <- DESeq(dds)

resultsNames(dds)
# [1] "Intercept"                 "treatment_OA_vs_AM"       
# [3] "treatment_OW_vs_AM"        "treatment_OWA_vs_AM"      
# [5] "generation_F4_vs_F0"       "treatmentOA.generationF4" 
# [7] "treatmentOW.generationF4"  "treatmentOWA.generationF4"

vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("treatment","generation"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=treatment,shape=generation)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

### TREATMENT EFFECT ###
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ generation + treatment)
dds <- DESeq(dds, test="LRT", reduced=~generation)
resultsNames(dds)
# [1] "Intercept"           "generation_F4_vs_F0" "treatment_OA_vs_AM" 
# [4] "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

resTreat <- results(dds, alpha = 0.05)
resTreat <- resTreat[order(resTreat$padj),] #ordering rows by adjusted p-value
head(resTreat)
# log2 fold change (MLE): treatment OWA vs AM 
# LRT p-value: '~ generation + treatment' vs '~ generation' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN1354_c0_g1    361.624      0.8228978  0.141554   99.9853
# TRINITY_DN73392_c0_g1   294.598     -0.7629778  0.112618   84.3020
# TRINITY_DN2266_c0_g1    569.405     -2.5441023  0.293738   82.0073
# TRINITY_DN1190_c0_g1    316.871      0.6053915  0.105221   73.8592
# TRINITY_DN34509_c0_g1   173.486      0.2121126  0.192976   71.0835
# TRINITY_DN9562_c0_g1    139.358     -0.0473908  0.141846   68.2358
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN1354_c0_g1  1.56555e-21 6.10157e-17
# TRINITY_DN73392_c0_g1 3.66412e-18 7.14027e-14
# TRINITY_DN2266_c0_g1  1.13869e-17 1.47931e-13
# TRINITY_DN1190_c0_g1  6.36200e-16 6.19881e-12
# TRINITY_DN34509_c0_g1 2.50166e-15 1.94999e-11
# TRINITY_DN9562_c0_g1  1.01852e-14 6.61597e-11
summary(resTreat)
# out of 128958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 886, 0.69%
# LFC < 0 (down)     : 643, 0.5%
# outliers [1]       : 12915, 10%
# low counts [2]     : 77069, 60%
# (mean count < 39)

resTreat <- resTreat[!is.na(resTreat$padj),] #filters data table to exclude NA results

degsTreat <- row.names(resTreat[resTreat$padj < 0.05,]) #differentially expressed genes from the environment
length(degsTreat)
#1529

## Generation effect
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ treatment + generation)
dds <- DESeq(dds, test="LRT", reduced=~treatment)
resultsNames(dds)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM" 
# [4] "treatment_OWA_vs_AM" "generation_F4_vs_F0"

resGen <- results(dds, alpha = 0.05)
resGen <- resGen[order(resGen$padj),] #ordering rows by adjusted p-value
head(resGen)
# log2 fold change (MLE): generation F4 vs F0 
# LRT p-value: '~ treatment + generation' vs '~ treatment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN8561_c0_g1    643.707      -0.892243 0.0805061  121.3638
# TRINITY_DN638_c0_g1     794.684      -2.435026 0.2120599  118.1561
# TRINITY_DN10760_c0_g1  5203.000       1.951769 0.1861498  102.2618
# TRINITY_DN34136_c0_g1  2191.023       2.004511 0.1951833   96.9525
# TRINITY_DN8105_c0_g1   2971.692      -2.002131 0.1962422   95.1547
# TRINITY_DN492_c1_g1    6580.513       1.320718 0.1394584   86.6035
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN8561_c0_g1  3.18105e-28 3.63642e-23
# TRINITY_DN638_c0_g1   1.60268e-27 9.16051e-23
# TRINITY_DN10760_c0_g1 4.86486e-24 1.85376e-19
# TRINITY_DN34136_c0_g1 7.10097e-23 2.02937e-18
# TRINITY_DN8105_c0_g1  1.76074e-22 4.02558e-18
# TRINITY_DN492_c1_g1   1.32614e-20 2.52663e-16
summary(resGen)
# out of 128958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 742, 0.58%
# LFC < 0 (down)     : 1931, 1.5%
# outliers [1]       : 12915, 10%
# low counts [2]     : 1728, 1.3%
# (mean count < 1)

resGen <- resGen[!is.na(resGen$padj),] #filters data table to exclude NA results

degsGen <- row.names(resGen[resGen$padj < 0.05,]) #differentially expressed genes from the environment
length(degsGen) 
#2673

#Interaction
## Generation effect
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ treatment + generation + treatment:generation)
dds <- DESeq(dds, test="LRT", reduced=~treatment + generation)
resultsNames(dds)
# [1] "Intercept"                 "treatment_OA_vs_AM"       
# [3] "treatment_OW_vs_AM"        "treatment_OWA_vs_AM"      
# [5] "generation_F4_vs_F0"       "treatmentOA.generationF4" 
# [7] "treatmentOW.generationF4"  "treatmentOWA.generationF4"

resInt<- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),] #ordering rows by adjusted p-value
head(resInt)
# log2 fold change (MLE): treatmentOWA.generationF4 
# LRT p-value: '~ treatment + generation + treatment:generation' vs '~ treatment + generation' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN1306_c0_g1   5226.925       -5.16109  0.449803   164.337
# TRINITY_DN2173_c0_g1   1051.624       -4.40042  0.426002   151.791
# TRINITY_DN11685_c0_g1   247.985       -4.46029  0.419968   117.967
# TRINITY_DN10411_c0_g1  1246.477       -2.54227  0.298087   114.689
# TRINITY_DN11890_c0_g1   237.432       -5.98824  0.683334   114.576
# TRINITY_DN21466_c0_g2  1243.539       -5.15112  0.559833   107.859
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN1306_c0_g1  2.12330e-35 1.35736e-30
# TRINITY_DN2173_c0_g1  1.08259e-32 3.46033e-28
# TRINITY_DN11685_c0_g1 2.11480e-25 4.50642e-21
# TRINITY_DN10411_c0_g1 1.07419e-24 1.45265e-20
# TRINITY_DN11890_c0_g1 1.13618e-24 1.45265e-20
# TRINITY_DN21466_c0_g2 3.16993e-23 3.37740e-19

summary(resInt)
# out of 128958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 644, 0.5%
# LFC < 0 (down)     : 754, 0.58%
# outliers [1]       : 2580, 2%
# low counts [2]     : 62451, 48%
# (mean count < 16)


resInt <- resInt[!is.na(resInt$padj),] #filters data table to exclude NA results

degsInt <- row.names(resInt[resInt$padj < 0.05,]) #differentially expressed genes from the environment
length(degsInt) 
#1398


## VENN DIAGRAM ##
# Total
length(degsGen)  # 2673
length(degsTreat)  # 1529
length(degsInt)  # 1398

# Intersections
length(intersect(degsGen,degsTreat))  # 373
length(intersect(degsGen,degsInt))  # 150
length(intersect(degsInt,degsTreat))  # 79

intEL <- intersect(degsGen,degsTreat)
length(intersect(degsInt,intEL)) # 30

# Number unique
2673-373-150-30 #2120
1529-373-79-30 #1047
1398-150-79-30 #1139

# unique, present in 2, present in 3
fit1 <- euler(c("Generation" = 2120, "Treatment" = 1047, "Interaction" = 1139, "Generation&Treatment" = 373, "Generation&Interaction" = 150, "Interaction&Treatment" = 79, "Generation&Treatment&Interaction" = 30))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#Heatmaps

topgenes <- head(rownames(resGen),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resTreat),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)
