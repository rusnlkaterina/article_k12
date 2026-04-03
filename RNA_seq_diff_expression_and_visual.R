# Differential expression analysis

# 1. Install and load packages
BiocManager::install("DESeq2") # Analyzing RNA-seq data

library(DESeq2)
library(ggplot2) # for plotting
library(dplyr)

library(AnnotationDbi) # for converting gene IDs

BiocManager::install("org.EcK12.eg.db")
library(org.EcK12.eg.db)

# for visualisation
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
library(pheatmap)
library(RColorBrewer)
# for volcano plots formation 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

install.packages("ggpubr")
library(ggpubr)

library(clusterProfiler) # for functional enrichment analysis, particularly for pathway analysis, 
# such as Gene Ontology (GO) and KEGG pathway analysis

install.packages("enrichR") # for functional analysis
library(enrichR)

# for volcano plots formation 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

# 2. Main analysis
# Import feature counts table (1)
countdata <- read.table("/Users/rusnlkaterina/Downloads/all_.counts", header = TRUE, skip = 0, row.names = 1)

# Make sure ID's are correct
head(countdata)

# Divide the info into 2 sections: Dam- vs WT (1.1), AlkB- vs WT (1.2), Dam- vs AlkB- (1.3)
## awk '{print $1,$5,$6,$7,$8,$9,$10}' all_.counts > all_WT_Dam.counts
## awk '{print $1,$2,$3,$4,$8,$9,$10}' all.counts > all_WT_AlkB.counts
## awk '{print $1,$2,$3,$4,$5,$6,$7}' all_.counts > all_Dam_AlkB.counts
countdata_Dam_WT <- read.table("/Users/rusnlkaterina/Downloads/all_WT_Dam.counts", header = TRUE, skip = 0, row.names = 1)
head(countdata_Dam_WT)

countdata_AlkB_WT <- read.table("/Users/rusnlkaterina/Downloads/all_WT_AlkB.counts", header = TRUE, skip = 0, row.names = 1)
head(countdata_AlkB_WT)

countdata_Dam_AlkB <- read.table("/Users/rusnlkaterina/Downloads/all_Dam_AlkB.counts", header = TRUE, skip = 0, row.names = 1)
head(countdata_Dam_AlkB)

# Import metadata file (2) / create manually
metadata <- read.table("/Users/rusnlkaterina/Downloads/metadata.txt", header = TRUE, row.names = 1)
head(metadata)

# Divide metadata into the two parts. Dam- vs WT (2.1), AlkB- vs WT (2.2), Dam- vs AlkB- (2.3)
metadata_Dam_WT <- read.table("/Users/rusnlkaterina/Downloads/metadata_Dam_WT.txt", header = TRUE, row.names = 1)
metadata_AlkB_WT <- read.table("/Users/rusnlkaterina/Downloads/metadata_AlkB_WT.txt", header = TRUE, row.names = 1)
metadata_Dam_AlkB <- read.table("/Users/rusnlkaterina/Downloads/metadata_Dam_AlkB.txt", header = TRUE, row.names = 1)

# Make sure ID's are correct
head(metadata_Dam_WT)
head(metadata_AlkB_WT)
head(metadata_Dam_AlkB)

# Make DESeq2 object from counts and metadata
ddsMat_Dam_WT <- DESeqDataSetFromMatrix(countData = countdata_Dam_WT,
                                        colData = metadata_Dam_WT,
                                        design = ~Group)
## design : The design of the comparisons to use. 
## use (~) before the name of the column variable to compare
ddsMat_AlkB_WT <- DESeqDataSetFromMatrix(countData = countdata_AlkB_WT,
                                         colData = metadata_AlkB_WT,
                                         design = ~Group)

ddsMat_Dam_AlkB <- DESeqDataSetFromMatrix(countData = countdata_Dam_AlkB,
                                          colData = metadata_Dam_AlkB,
                                          design = ~Group)


# Find differential expressed genes
# Run DESEq2
ddsMat_Dam_WT <- DESeq(ddsMat_Dam_WT)
ddsMat_AlkB_WT <- DESeq(ddsMat_AlkB_WT)
ddsMat_Dam_AlkB <- DESeq(ddsMat_Dam_AlkB)

# Get basic statisics about the number of significant genes
res_Dam_WT <- results(ddsMat_Dam_WT, contrast = c("Group", "group1", "group2"), alpha = 0.05) # Dam- vs WT
res_AlkB_WT <- results(ddsMat_AlkB_WT, contrast = c("Group", "group1", "group3"), alpha = 0.05) # AlkB- vs WT
res_Dam_AlkB <- results(ddsMat_Dam_AlkB, contrast = c("Group", "group2", "group3"), alpha = 0.05) # Dam- vs AlkB-

# Generate summary of testing
summary(res_Dam_WT)
summary(res_AlkB_WT)
summary(res_Dam_AlkB)

# Subset for only significant genes (p < 0.05)
res_Dam_WT_sig <- subset(res_Dam_WT, padj < 0.05)
res_AlkB_WT_sig  <- subset(res_AlkB_WT, padj < 0.05)
res_Dam_AlkB_sig <- subset(res_Dam_AlkB, padj < 0.05)

head(res_Dam_WT_sig)
head(res_AlkB_WT_sig)
head(res_Dam_AlkB_sig)

nrow(res_Dam_WT_sig) # 2780 top dif. expressed genes # Dam- vs WT
nrow(res_AlkB_WT_sig) # 1981 top dif. expressed genes # AlkB- vs WT
nrow(res_Dam_AlkB_sig) # 2237 top dif. expressed genes # AlkB- vs WT

nrow(res_Dam_WT) # 4494 all the genes
nrow(res_AlkB_WT) # 4494 all the genes
nrow(res_Dam_AlkB) # 4494 all the genes

# 3. Visualisation
# Heatmap plot
# Dam_WT
# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_Dam_WT)$Group),
  Replicate = factor(colData(ddsMat_Dam_WT)$Replicate),
  row.names = colData(ddsMat_Dam_WT)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c("group1" = "pink", "group2" = "lightblue"),
  Replicate = c(Rep1 = "white", Rep2 = "grey", Rep3 ="black")
)

mat_50_Dam_AlkB <- assay(ddsMat_Dam_AlkB[rownames(tail(res_Dam_AlkB_sig[order(abs(res_Dam_AlkB_sig$log2FoldChange)),], 50))])
nrow(mat_50_Dam_AlkB)

mat_50_Dam_AlkB <- subset(mat_50_Dam_AlkB, is.na(row.names(mat_50_Dam_AlkB)) == FALSE) # delete NAs
nrow(mat_50_Dam_AlkB) # 50 - good

mat_50_Dam_WT <- assay(ddsMat_Dam_WT[rownames(tail(res_Dam_WT_sig[order(abs(res_Dam_WT_sig$log2FoldChange)),], 50))])
nrow(mat_50_Dam_WT)

mat_50_Dam_WT <- subset(mat_50_Dam_WT, is.na(row.names(mat_50_Dam_WT)) == FALSE) # delete NAs
nrow(mat_50_Dam_WT) # 50 - good

# Make heatmap with pheatmap function.
# see more in documentation for customization
pheatmap(mat = mat_50_Dam_WT, 
         color = colorRampPalette(brewer.pal(9, "Blues"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         fontsize = 8,
         show_colnames = F)

# save image 500 width x 850 height

# Dam_AlkB
annotation_col = data.frame(
  Group = factor(colData(ddsMat_Dam_AlkB)$Group),
  Replicate = factor(colData(ddsMat_Dam_AlkB)$Replicate),
  row.names = colData(ddsMat_Dam_AlkB)$sampleid
)

ann_colors = list(
  Group = c("group3" = "lightgreen", "group2" = "lightblue"),
  Replicate = c(Rep1 = "white", Rep2 = "grey", Rep3 ="black")
)


mat_50_Dam_AlkB <- assay(ddsMat_Dam_AlkB[rownames(tail(res_Dam_AlkB_sig[order(abs(res_Dam_AlkB_sig$log2FoldChange)),], 50))])
nrow(mat_50_Dam_AlkB)

mat_50_Dam_AlkB <- subset(mat_50_Dam_AlkB, is.na(row.names(mat_50_Dam_AlkB)) == FALSE) # delete NAs
nrow(mat_50_Dam_AlkB) # 50 - good

pheatmap(mat = mat_50_Dam_AlkB, 
         color = colorRampPalette(brewer.pal(9, "Oranges"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         fontsize = 8,
         show_colnames = F)

# AlkB_WT
annotation_col = data.frame(
  Group = factor(colData(ddsMat_AlkB_WT)$Group),
  Replicate = factor(colData(ddsMat_AlkB_WT)$Replicate),
  row.names = colData(ddsMat_AlkB_WT)$sampleid
)

ann_colors = list(
  Group = c("group1" = "pink", "group3" = "lightgreen"),
  Replicate = c(Rep1 = "white", Rep2 = "grey", Rep3 ="black")
)

mat_50_AlkB_WT <- assay(ddsMat_AlkB_WT[rownames(tail(res_AlkB_WT_sig[order(abs(res_AlkB_WT_sig$log2FoldChange)),], 50))])
nrow(mat_50_AlkB_WT)

mat_50_AlkB_WT <- subset(mat_50_AlkB_WT, is.na(row.names(mat_50_AlkB_WT)) == FALSE) # delete NAs
nrow(mat_50_AlkB_WT) # 50 - good

pheatmap(mat = mat_50_AlkB_WT, 
         color = colorRampPalette(brewer.pal(9, "Greens"))(255), 
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors, 
         fontsize = 8,
         show_colnames = F)

# Create a volcano plot
write.csv(res_Dam_WT_sig, "/Users/rusnlkaterina/Downloads/res_Dam_WT_sig.csv", row.names = TRUE)
write.csv(res_AlkB_WT_sig, "/Users/rusnlkaterina/Downloads/res_AlkB_WT_sig.csv", row.names = TRUE)
write.csv(res_Dam_AlkB_sig, "/Users/rusnlkaterina/Downloads/res_Dam_AlkB_sig.csv", row.names = TRUE)

# Import DGE (differential gene expression) results
df_Dam_WT <- read.csv("/Users/rusnlkaterina/Downloads/res_Dam_WT_sig.csv", header = TRUE)
head(df_Dam_WT)

df_AlkB_WT <- read.csv("/Users/rusnlkaterina/Downloads/res_AlkB_WT_sig.csv", header = TRUE)
head(df_AlkB_WT)

df_Dam_AlkB <- read.csv("/Users/rusnlkaterina/Downloads/res_Dam_AlkB_sig.csv", header = TRUE)
head(df_Dam_AlkB)

selected_genes <- c('topA', 'topB', 'rnhA', 'rnhB', 'gyrA', 'gyrB', 'dinG', 'recQ', 'recG',
                    'recD', 'recB', 'recA', 'recF', 'uvsW', 'pcrA', 'uvrB', 'radD', 'rho',
                    'nus6', 'nusA', 'cas3', 'uvrD', 'priA')

# Dam_WT
df_Dam_WT$log2FoldChange <- df_Dam_WT$log2FoldChange * -1 # to make it appropriate for y-axis on the plot

# create custom key-value pairs for 'up', 'down' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  df_Dam_WT$log2FoldChange < 0, 'red2',
  ifelse(df_Dam_WT$log2FoldChange > 0, 'royalblue',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red2'] <- 'down'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'up'

EnhancedVolcano(df_Dam_WT,
                lab = df_Dam_WT$X,
                x = 'log2FoldChange',
                xlim = c(-3,3),
                y = 'padj',
                selectLab = selected_genes,
                colAlpha = 1/5,
                pointSize = 1.5,
                labSize = 4.0,
                boxedLabels = TRUE,
                drawConnectors=TRUE,
                colCustom = keyvals,
                pCutoff = 0.05, # Adjusted p-value cutoff
                FCcutoff = 0.0) # Log2 fold change cutoff

# AlkB_WT
df_AlkB_WT$log2FoldChange <- df_AlkB_WT$log2FoldChange * -1

# create custom key-value pairs for 'up', 'down' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  df_AlkB_WT$log2FoldChange < 0, 'red2',
  ifelse(df_AlkB_WT$log2FoldChange > 0, 'green3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red2'] <- 'down'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'green3'] <- 'up'

EnhancedVolcano(df_AlkB_WT,
                lab = df_AlkB_WT$X,
                x = 'log2FoldChange',
                xlim = c(-3,3),
                y = 'padj',
                selectLab = selected_genes,
                colAlpha = 1/5,
                pointSize = 1.5,
                labSize = 4.0,
                boxedLabels = TRUE,
                drawConnectors=TRUE,
                colCustom = keyvals,
                pCutoff = 0.05, # Adjusted p-value cutoff
                FCcutoff = 0.0) # Log2 fold change cutoff

# Dam_AlkB
df_Dam_AlkB$log2FoldChange <- df_Dam_AlkB$log2FoldChange * -1

keyvals <- ifelse(
  df_Dam_AlkB$log2FoldChange < 0, 'royalblue',
  ifelse(df_Dam_AlkB$log2FoldChange > 0, 'green3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'royalblue'] <- 'down'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'green3'] <- 'up'

EnhancedVolcano(df_Dam_AlkB,
                lab = df_Dam_AlkB$X,
                x = 'log2FoldChange',
                xlim = c(-5,5),
                y = 'padj',
                selectLab = selected_genes,
                colAlpha = 1/5,
                pointSize = 1.5,
                labSize = 3.0,
                boxedLabels = TRUE,
                drawConnectors=TRUE,
                colCustom = keyvals,
                pCutoff = 0.05, # Adjusted p-value cutoff
                FCcutoff = 0.0) # Log2 fold change cutoff

# Boxplots
file <- read.csv("/Users/rusnlkaterina/Downloads/boxplot_data_wt.csv", header = TRUE)
order <- c("Non R-loops containing", "R-loops containing", "m6A/6mA and R-loops containing") 

# for statistical significance test
my_comparisons <- list(
  c("Non R-loops containing", "R-loops containing"),
  c("R-loops containing", "m6A/6mA and R-loops containing"),
  c("Non R-loops containing", "m6A/6mA and R-loops containing")
)

file %>% 
  ggplot(aes(x=factor(subtype, level=order), y=rpkm)) +
  geom_boxplot(outlier.shape = NA) + # to not display outliers (since there is jitter / as they will be shown by geom_jitter)
  geom_jitter(width=0.34, alpha=0.07) +
  geom_point(colour = 'salmon', size = 0.5)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = 4.5, method = "kruskal.test") + # kruskal.test better for wilcox.test
  scale_y_log10() +
  theme_minimal()