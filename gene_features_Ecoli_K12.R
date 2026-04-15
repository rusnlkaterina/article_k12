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

txdb <- makeTxDbFromGFF(file="/Users/rusnlkaterina/Downloads/ecoli-k12_strain.gtf")
txdb

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-300,300))
peakAnnoList
plotAnnoBar(peakAnnoList)

# Boxplots formation (RNA-seq and DRIP-seq intersection analysis)
# Gene names retrieval 
S9.6_WT_annot <- data.frame(peakAnnoList[["S9.6 WT"]]@anno)
nrow(S9.6_WT_annot)

S9.6_WT_gene_ids <- S9.6_WT_annot$geneId

write.csv(S9.6_WT_gene_ids, "/Users/rusnlkaterina/Downloads/S9.6_WT_gene_ids.csv", row.names=FALSE)
# remove first line with "x"
# grep -f S9.6_WT_gene_ids.csv ecoli-k12_strain.gtf > extracted_rows.txt
# awk -F'gene "|"; locus_tag' '{print $2}' extracted_rows.txt > extracted_gene_names.txt
# grep -v 'gene_biotype' extracted_gene_names.txt > extracted_gene_namess.txt
# sort extracted_gene_namess.txt | uniq > extracted_gene_namess_unique.txt