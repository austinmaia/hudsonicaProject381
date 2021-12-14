## Import or install the libraries that we're likely to need
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 
library(eulerr)
setwd("~/Documents/GitHub/hudsonicaProject381/data")

# Import the counts matrix
countsTable <- read.table("salmon.gene.counts.f0.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
#[1] 130580    12
countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("F0.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

### DESeq
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)

# Filter out genes with too few reads 
dds <- dds[rowSums(counts(dds)) >120]
dim(dds)
# [1] 81694    12

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)
resultsNames(dds)

### PCA to visualize
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd, intgroup="treatment", returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2, color=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Comparison F0: AM vs OW

resOW <- results(dds, contrast=c("treatment","AM","OW"), alpha=0.05)

resOW <- resOW[order(resOW$padj),]
head(resOW)
summary(resOW)


# Comparison AM vs OA


resOA <- results(dds, contrast=c("treatment","AM","OA"), alpha=0.05)
resOA <- resOA[order(resOA$padj),]
head(resOA)
summary(resOA)


# Comparison AM vs OWA


resOWA <- results(dds, contrast=c("treatment","AM","OWA"), alpha=0.05)
resOWA <- resOWA[order(resOWA$padj),]
head(resOWA)
summary(resOWA)

### DEGS and Venn ### 

names(resOWA)[6] <- 'pOWA'
names(resOA)[6] <- 'pOA'
names(resOW)[6] <- 'pOW'

fixed <- merge.data.frame(resOWA, resOA, by=0, all=TRUE)
fixed$pOW <- resOW$pOW
row.names(fixed) <- fixed$Row.names
fixed <- fixed[!is.na(fixed$pOA),]
fixed <- fixed[!is.na(fixed$pOW),]
fixed <- fixed[!is.na(fixed$pOWA),]
degsCombo <- row.names(fixed[fixed$pOWA < 0.05,])
length(degsCombo) #4243

degsOA <- row.names(fixed[fixed$pOA < 0.05,])
length(degsOA) #1357

degsOW <- row.names(fixed[fixed$pOW < 0.05,])
length(degsOW) #7975

## VENN DIAGRAM ##
# Total
length(degsOA)  # 1357
length(degsOW)  # 7975
length(degsCombo)  # 4243

# Intersections
length(intersect(degsOA,degsOW))  # 150
length(intersect(degsOA,degsCombo))  # 575
length(intersect(degsCombo,degsOW))  # 479

intEL <- intersect(degsOA,degsOW)
length(intersect(degsCombo,intEL)) # 70

# Number unique
1357-150-575-70 #OA 562
7975-150-479-70 #OW 7276
4243-575-479-70 #Combo 3119

# unique, present in 2, present in 3
fit1 <- euler(c("OA" = 562, "OW" = 7276, "Combined" = 3119, "OA&OW" = 150, "OA&Combined" = 575, "OW&Combined" = 479, "OA&OW&Combined" = 70))

plot(fit1,  lty = 1:3, quantities = TRUE)


## Heat map
library(pheatmap)

topgenes <- head(rownames(resOA),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resOWA),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resOW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)
