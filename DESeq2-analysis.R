# SCRIPT BY DR. BRIALLEN LOBB

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(stringr)
library(GenomicFeatures)
library(reshape2)
library(viridis)

#H_sapiens_A549
#DESeq2
samples <- list.files(path = "~/Desktop/Salmon_quants/H_sapiens/", full.names = T)
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "/Users/briallenlobb/Desktop/Salmon_quants/H_sapiens//", "") %>% str_replace("_transcripts_quant", "")
files <- files[c(1:3,5:7,9:16)]
txdb <- makeTxDbFromGFF( file = "~/Desktop/GCF_000001405.40_H_sapiens_GRCh38.p14_genomic.gff")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select( txdb, k, "GENEID", "TXNAME" )
txi <- tximport(files[1:6], type="salmon", countsFromAbundance="lengthScaledTPM",tx2gene = tx2gene)
samples <- as.data.frame(names(files)[1:6])
samples$condition <- factor(c(rep("Not",3),(rep("IFN",3))))
rownames(samples) <- samples[,1]
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Not")
hA549dds <- DESeq(dds)
hA549res <- results(hA549dds)
write.table(hA549res,"~/Desktop/H_sapiens_A549_DESeq2_results.txt",quote = FALSE,sep="\t")
#PCA
hA549vsd <- vst(hA549dds, blind=FALSE)
pcaData <- plotPCA(hA549vsd, intgroup=c("condition"), returnData=TRUE)
pcaData$name <- sapply(strsplit(pcaData$name,"_"),"[",1)
rownames(pcaData) <- pcaData$name
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- samples$condition
pdf("~/Desktop/H_sapiens_A549_PCA.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() + geom_text_repel()
dev.off()
#EXPORT RESULTS
hA549res05 <- results(hA549dds, alpha=0.05)
summary(hA549res05)
write.table(hA549res05[!is.na(hA549res05$padj) & hA549res05$padj < 0.05 & hA549res05$log2FoldChange > 0,], "~/Desktop/H_sapiens_A549_DEG_up.tsv", quote = FALSE,sep="\t")
write.table(hA549res05[!is.na(hA549res05$padj) & hA549res05$padj < 0.05 & hA549res05$log2FoldChange < 0,], "~/Desktop/H_sapiens_A549_DEG_down.tsv", quote = FALSE,sep="\t")

#H_sapiens_RPTEC
samples <- list.files(path = "~/Desktop/Salmon_quants/H_sapiens/", full.names = T)
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "/Users/briallenlobb/Desktop/Salmon_quants/H_sapiens//", "") %>% str_replace("_transcripts_quant", "")
files <- files[c(9:10,12:14,16)]
txi <- tximport(files, type="salmon", countsFromAbundance="lengthScaledTPM",tx2gene = tx2gene)
samples <- as.data.frame(names(files))
samples$condition <- factor(c(rep("Not",3),(rep("IFN",3))))
rownames(samples) <- samples[,1]
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Not")
hRPTECdds <- DESeq(dds)
hRPTECres <- results(hRPTECdds)
write.table(hRPTECres,"~/Desktop/H_sapiens_RPTEC_DESeq2_results.txt",quote = FALSE,sep="\t")
#PCA
hRPTECvsd <- vst(hRPTECdds, blind=FALSE)
pcaData <- plotPCA(hRPTECvsd, intgroup=c("condition"), returnData=TRUE)
pcaData$name <- sapply(strsplit(pcaData$name,"_"),"[",1)
rownames(pcaData) <- pcaData$name
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- samples$condition
pdf("~/Desktop/H_sapiens_RPTEC_PCA.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() + geom_text_repel()
dev.off()
#EXPORT RESULTS
hRPTECres05 <- results(hRPTECdds, alpha=0.05)
summary(hRPTECres05)
write.table(hRPTECres05[!is.na(hRPTECres05$padj) & hRPTECres05$padj < 0.05 & hRPTECres05$log2FoldChange > 0,], "~/Desktop/H_sapiens_RPTEC_DEG_up.tsv", quote = FALSE,sep="\t")
write.table(hRPTECres05[!is.na(hRPTECres05$padj) & hRPTECres05$padj < 0.05 & hRPTECres05$log2FoldChange < 0,], "~/Desktop/H_sapiens_RPTEC_DEG_down.tsv", quote = FALSE,sep="\t")

#E_fuscus
#DESeq2
samples <- list.files(path = "~/Desktop/Salmon_quants/E_fuscus/", full.names = T)
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "/Users/briallenlobb/Desktop/Salmon_quants/E_fuscus//", "") %>% str_replace("_transcripts_quant", "")
txdb <- makeTxDbFromGFF( file = "~/Desktop/GCF_027574615.1_E_fuscus_DD_ASM_mEF_20220401_genomic.gff")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select( txdb, k, "GENEID", "TXNAME" )
txi <- tximport(files, type="salmon", countsFromAbundance="lengthScaledTPM",tx2gene = tx2gene)
samples <- as.data.frame(names(files))
samples$condition <- factor(c(rep("Not",3),(rep("IFN",3))))
rownames(samples) <- samples[,1]
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Not")
efdds <- DESeq(dds)
efres <- results(efdds)
write.table(efres,"~/Desktop/E_fuscus_DESeq2_results.txt",quote = FALSE,sep="\t")
#PCA
efvsd <- vst(efdds, blind=FALSE)
pcaData <- plotPCA(efvsd, intgroup=c("condition"), returnData=TRUE)
pcaData$name <- sapply(strsplit(pcaData$name,"_"),"[",1)
rownames(pcaData) <- pcaData$name
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- samples$condition
pdf("~/Desktop/E_fuscus_PCA.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() + geom_text_repel()
dev.off()
#EXPORT RESULTS
efres05 <- results(efdds, alpha=0.05)
summary(efres05)
write.table(efres05[!is.na(efres05$padj) & efres05$padj < 0.05 & efres05$log2FoldChange > 0,], "~/Desktop/E_fuscus_DEG_up.tsv", quote = FALSE,sep="\t")
write.table(efres05[!is.na(efres05$padj) & efres05$padj < 0.05 & efres05$log2FoldChange < 0,], "~/Desktop/E_fuscus_DEG_down.tsv", quote = FALSE,sep="\t")

#P_alecto
#DESeq2
samples <- list.files(path = "~/Desktop/Salmon_quants/P_alecto/", full.names = T)
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "/Users/briallenlobb/Desktop/Salmon_quants/P_alecto//", "") %>% str_replace("_transcripts_quant", "")
txdb <- makeTxDbFromGFF( file = "~/Desktop/GCF_000325575.1_P_alecto_ASM32557v1_genomic.gff")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select( txdb, k, "GENEID", "TXNAME" )
files <- files[c(1:3,5:7)]
txi <- tximport(files, type="salmon", countsFromAbundance="lengthScaledTPM",tx2gene = tx2gene)
samples <- as.data.frame(names(files))
samples$condition <- factor(c(rep("Not",3),(rep("IFN",3))))
rownames(samples) <- samples[,1]
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Not")
padds <- DESeq(dds)
pares <- results(padds)
write.table(pares,"~/Desktop/P_alecto_DESeq2_results.txt",quote = FALSE,sep="\t")
#PCA
pavsd <- vst(padds, blind=FALSE)
pcaData <- plotPCA(pavsd, intgroup=c("condition"), returnData=TRUE)
pcaData$name <- sapply(strsplit(pcaData$name,"_"),"[",1)
rownames(pcaData) <- pcaData$name
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- samples$condition
pdf("~/Desktop/P_alecto_PCA.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() + geom_text_repel()
dev.off()
#EXPORT RESULTS
pares05 <- results(padds, alpha=0.05)
summary(pares05)
write.table(pares05[!is.na(pares05$padj) & pares05$padj < 0.05 & pares05$log2FoldChange > 0,], "~/Desktop/P_alecto_DEG_up.tsv", quote = FALSE,sep="\t")
write.table(pares05[!is.na(pares05$padj) & pares05$padj < 0.05 & pares05$log2FoldChange < 0,], "~/Desktop/P_alecto_DEG_down.tsv", quote = FALSE,sep="\t")
