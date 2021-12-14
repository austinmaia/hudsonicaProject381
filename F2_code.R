#F2 data
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 
library(pheatmap)

setwd("~/Documents/GitHub/hudsonicaProject381/data")

countsTable <- read.table("salmon.gene.counts.f2.txt", header=TRUE, row.names = 1)
head(countsTable)
countsTableRound <-round(countsTable)
conds <- read.delim("F2.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ treatment)
dim(dds)
dds <- dds[rowSums(counts(dds)) > 80]

dds <- DESeq(dds)
vsd <- vst(dds, blind=F)


#list results generated
resultsNames(dds)
#[1] "Intercept"          "treatment_OA_vs_AM" "treatment_OW_vs_AM"
resF2OA_OW <- results(dds, contrast=c("treatment","OA","OW"), alpha=0.05)
resF2OA_OW <- resF2OA_OW[order(resF2OA_OW$padj),]
head(resF2OA_OW)
summary(resF2OA_OW)
resF2OA_OW <- resF2OA_OW[!is.na(resF2OA_OW$padj),]
degsF2OA_OW <- row.names(resF2OA_OW[resF2OA_OW$padj < 0.05,])

resF2OW <- results(dds,contrast=c("treatment","OW","AM"), alpha=0.05)
resF2OW <- resF2OW[order(resF2OW$padj),]
head(resF2OW)
summary(resF2OW)
resF2OW <- resF2OW[!is.na(resF2OW$padj),]
degsF2OW <- row.names(resF2OW[resF2OW$padj < 0.05,])

resF2OA <- results(dds, contrast=c("treatment","OA","AM"), alpha=0.05)
resF2OA <- resF2OA[order(resF2OA$padj),]
head(resF2OA)
summary(resF2OA)
resF2OA <- resF2OA[!is.na(resF2OA$padj),]
degsF2OA <- row.names(resF2OA[resF2OA$padj < 0.05,])

#looking at specific gene compared across treatment
plotCounts(dds, gene="TRINITY_DN45760_c0_g1", intgroup="treatment")

#making venni diagram
library(eulerr)

# Total
length(degsF2OA)  # 379
length(degsF2OW)  # 1472
length(degsF2OA_OW)  #1912 

# Intersections
length(intersect(degsF2OA,degsF2OW))  # 195
length(intersect(degsF2OA,degsF2OA_OW))  # 137
length(intersect(degsF2OW,degsF2OA_OW))  # 832

intEL <- intersect(degsF2OA,degsF2OW)
length(intersect(degsF2OA_OW,intEL)) # 2
# Number unique
379-137-195-2 #45
1472-195-832-2 #443
1912-137-832-2 #941

fit1 <- euler(c("OA"=45, "OW"=443, "OA_OW"=941, "OA&OW"=195, "OW&OA_OW"=137, "OA&OA_OW"=832, "OA&OW&OA_OW"=2))

plot(fit1,  lty = 1:3, quantities = TRUE)
              
plot(fit1, quantities = TRUE, fill = "transparent",lty = 1:3,labels = list(font = 4))
   

#edits for consistency#
names(resF2OA_OW)[6] <- 'pOWA'
names(resF2OA)[6] <- 'pOA'
names(resF2OW)[6] <- 'pOW'

fixed <- merge.data.frame(resF2OA_OW, resF2OA, by=0, all=TRUE)
fixed$pOW <- resF2OW$pOW
row.names(fixed) <- fixed$Row.names
fixed <- fixed[!is.na(fixed$pOA),]
fixed <- fixed[!is.na(fixed$pOW),]
fixed <- fixed[!is.na(fixed$pOWA),]
degsCombo <- row.names(fixed[fixed$pOWA < 0.05,])
length(degsCombo) #2167

degsOA <- row.names(fixed[fixed$pOA < 0.05,])
length(degsOA) #406

degsOW <- row.names(fixed[fixed$pOW < 0.05,])
length(degsOW) #1708

## VENN DIAGRAM ##
# Total
length(degsOA)  # 406
length(degsOW)  # 1708
length(degsCombo)  # 2167

# Intersections
length(intersect(degsOA,degsOW))  # 10
length(intersect(degsOA,degsCombo))  # 148
length(intersect(degsCombo,degsOW))  # 69

intEL <- intersect(degsOA,degsOW)
length(intersect(degsCombo,intEL)) # 2

# Number unique
406-10-148-2 #OA 246
1708-10-69-2 #OW 1627
2167-148-69-2 #Combo 1948

# unique, present in 2, present in 3
fit1 <- euler(c("OA" = 246, "OW" = 1627, "OA vs OW" = 1948, "OA&OW" = 10, "OA&OA vs OW" = 148, "OW&OA vs OW" = 69, "OA&OW&OA vs OW" = 2))

plot(fit1,  lty = 1:3, quantities = TRUE)

fit2 <- euler(c("OA" = 246, "OW" = 1627, "OA&OW" = 10))

plot(fit2,  lty = 1:3, quantities = TRUE)

#Heatmaps
#Heatmaps

topgenes <- head(rownames(resF2OA),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF2OA_OW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF2OW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)
