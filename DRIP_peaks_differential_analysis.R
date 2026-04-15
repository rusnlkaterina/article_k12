# Install packages
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DiffBind")

# Load libraries
library(DiffBind)
library(tidyverse)
library(rtracklayer)

## S9.6 WT vs S9.6 Dam-
# Create a sample data frame
my_data <- data.frame(
  SampleID = c("WT_S9.6_1", "WT_S9.6_2", "Dam_min_S9.6_1", "Dam_min_S9.6_2"),
  Condition = c("WT", "WT", "KO", "KO"), # KO means knockout
  Replicate = c("1","2","1","2"),
  bamReads = c("/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_S96_1_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_S96_2_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_S9_6_1_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_S9_6_2_aln.bam"),
  ControlID = c("WT", "WT", "Dam_min", "Dam_min"),
  bamControl = c("/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_E_Coli_aln.bam"),
  Peaks = c("/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_S96.bed",
            "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/wt_S96.bed",
            "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_S9_6.bed",
            "/Users/rusnlkaterina/Downloads/WT_Dam_min_S96/DAM_min_S9_6.bed"),
  PeakCaller = c("bed", "bed", "bed", "bed")
)

## S9.6 WT vs S9.6 AlkB-
# Create a sample data frame
my_data <- data.frame(
  SampleID = c("WT_S9.6_1", "WT_S9.6_2", "AlkB_min_S9.6_1", "AlkB_min_S9.6_2"),
  Condition = c("WT", "WT", "KO", "KO"), # KO means knockout
  Replicate = c("1","2","1","2"),
  bamReads = c("/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_S96_1_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_S96_2_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_S9_6_1_aln.bam",
               "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_S9_6_2_aln.bam"),
  ControlID = c("WT", "WT", "AlkB_min", "AlkB_min"),
  bamControl = c("/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_E_Coli_aln.bam",
                 "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_E_Coli_aln.bam"),
  Peaks = c("/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_S96.bed",
            "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/wt_S96.bed",
            "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_S9_6.bed",
            "/Users/rusnlkaterina/Downloads/WT_AlkB_min_S96/AlkB_min_S9_6.bed"),
  PeakCaller = c("bed", "bed", "bed", "bed")
)

# Export the data frame to a CSV file
# The file will be saved in your current working directory unless a full path is specified.
write.csv(my_data, "/Users/rusnlkaterina/Downloads/metadata_peaks_wt_alkb.csv", row.names = FALSE) # or metadata_peaks_wt_dam.csv


# Read a csv file
samples <- read.csv("/Users/rusnlkaterina/Downloads/metadata_peaks_wt_alkb.csv")

# Look at the loaded metadata
names(samples)          
samples

metadata_peaks <- dba(sampleSheet="/Users/rusnlkaterina/Downloads/metadata_peaks_wt_alkb.csv")
metadata_peaks

metadata_peaks.counted <- dba.count(metadata_peaks, summits=500)
metadata_peaks.counted <- dba.contrast(metadata_peaks.counted, categories=DBA_CONDITION,
                                       minMembers=2)

# DESeq2
metadata_peaks.analysed <- dba.analyze(metadata_peaks.counted)

dba.show(metadata_peaks.analysed, bContrasts=T)

dba.plotVolcano(metadata_peaks.analysed, bUsePval = TRUE, th=0.05)

report <- dba.report(metadata_peaks.analysed)
report

## S9.6 WT vs S9.6 Dam- / S9.6 WT vs S9.6 AlkB-
# Rename a single column
colnames(mcols(report))[colnames(mcols(report)) == "p-value"] <- "pvalue"
report

with(report, plot(Fold, -log10(pvalue), pch=20, main="Volcano plot"))

write.csv(report, "/Users/rusnlkaterina/Downloads/inf_wt_dam.csv", row.names=FALSE) # or inf_wt_alkb.csv

# Load libraries [2]
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(GenomicFeatures)

# Peaks annotation
# Load data
samplefiles <- list.files("/Users/rusnlkaterina/Downloads/inf_wt_dam", # or /inf_wt_alkb p.s. before convert .csv into .tsv manually
                          pattern= ".tsv", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("inff") 

txdb <- makeTxDbFromGFF(file="/Users/rusnlkaterina/Downloads/ecoli-k12_strain.gtf")
txdb

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-300,300))
peakAnnoList
plotAnnoBar(peakAnnoList)

# Gene names retrieval 
inf <- data.frame(peakAnnoList[["inff"]]@anno)
nrow(inf)

inf_gene_ids <- inf$geneId

write.csv(inf_gene_ids, "/Users/rusnlkaterina/Downloads/inf_gene_ids_WT_Dam.csv", row.names=FALSE) # or inf_gene_ids_WT_AlkB.csv
# sed 's/"//g' inf_gene_ids_WT_Dam.csv > inf_gene_ids_WT_Dam_.csv
# Then, this resulted file should be implemented in the python script 'RNA-seq_K12.ipynb'.
file_path <- "/Users/rusnlkaterina/Downloads/inf_color_final_wt_alkb.txt" # OR inf_color_final_wt_dam.txt
words_vector <- scan(file_path, what = character(), sep = " ")
report$color <- words_vector # for WT / AlkB- OR WR / Dam-

# Plot with colors based on the 'group' column
plot(report$Fold, -log10(report$pvalue), col = report$color, pch = 19, main = "Scatter Plot with Colored Dots")

# count number of dots:
# grep -o -w "black" inf_color_final_dam_alkb.txt | wc -l # OR green OR red OR blue (it depends on the dataset)