## DIP/DRIP-seq analysis of E.coli
# gene features annotation

BiocManager::install("ChIPseeker")

# Load libraries
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")

library(GenomicFeatures)

# Broad peaks annotation

# Load data
samplefiles <- list.files("/Users/rusnlkaterina/Downloads/peaks_ecoli_k12_wt_dam",
                          pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("m6A/6mA Dam-", "S9.6 Dam-",
                        "m6A/6mA WT", "S9.6 WT") 

txdb <- makeTxDbFromGFF(file="/Users/rusnlkaterina/Downloads/ecoli-k12_strain.gff3")
txdb

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-300,300))
peakAnnoList
plotAnnoBar(peakAnnoList)