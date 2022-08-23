library(DESeq2)
library(biomaRt)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

setwd("path/to/wd")
list.files("read_counts/")

# Load the sample data from csv file
coldata <- read.csv("docs/sample_layout.csv", header=T)
coldata

# Creat raw counts matrix
count_files <- list.files("read_counts/")[grep("counts$", list.files("read_counts/"))]
sample_names <- gsub("_R1.trim.sorted.counts.counts", "", count_files)

file_index <- match(coldata$FileID, sample_names)
count_files <- count_files[file_index]
dge <- readDGE(files = count_files,
               columns = c(1, 7),
               path = "read_counts",
               skip = 1)
cts <- dge$counts
colnames(cts) <- gsub("_R1.trim.sorted.counts", "", 
                      colnames(cts))
colnames(cts)

# Create DESeq2 data set object from the count matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)
# Save raw unfiltered counts
write.csv(counts(dds), file = "counts_raw_unfiltered.csv", 
          row.names = T)
save(dds, file = "DESeqOBJ.RData")
load("expression_data/DESeqOBJ.RData")

# Get annotation 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
listDatasets(ensembl)$dataset[grep("aries",listDatasets(ensembl)$dataset)]

searchDatasets(mart = ensembl, pattern = "aries")
ensembl <- useDataset("oaries_gene_ensembl", mart=ensembl)
# Check for attribute names
listAttributes(ensembl)$name[grep("symbol",listAttributes(ensembl)$name)]
listAttributes(ensembl)$name[grep("entrez",listAttributes(ensembl)$name)]
listAttributes(ensembl)$name[grep("type",listAttributes(ensembl)$name)]

attributes <- c("ensembl_gene_id", "uniprot_gn_symbol", 
                "entrezgene_id", "description", "gene_biotype")

ids <- rownames(counts(dds))
annot <- getBM(attributes = attributes, values = ids, mart = ensembl)
annot <- annot[!duplicated(annot$ensembl_gene_id),]
write.csv(annot, file = "gene_annotation.csv")

###############################################################
# Sample clustering and visualization
vsd <- vst(dds)
head(assay(vsd))

df <- as.data.frame(colData(dds)[,c("Organ","Group")])

# Get top 500 genes with the highest variance
vars <- rowVars(assay(vsd))
names(vars) <- rownames(assay(vsd))
vars <- sort(vars, decreasing=TRUE)
vars <- vars[1:500]
topVars <- assay(vsd)[names(vars),]

colors <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                            "RdYlBu")))(100)
tiff("all_samples_heatmap_euclidean_wardD2.tiff", width=600, 
     height = 600)
pheatmap(topVars, scale = "row", fontsize_row=3,
         color = colors,
         annotation_col=df, 
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         labels_col = coldata$SampleID,
         show_rownames = F)
dev.off()

## Plot heatmap of sample-to-sample distances based on 
## top 500 genes by variance
sampleDists <- dist(t(topVars), method = "euclidean")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$Organ, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff("sample_distance_matrix_complete.tiff", width = 600, height = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         fontsize_row = 5,
         clustering_method = "complete")
dev.off()

## Create PCA plot for all samples
pca_dat <- prcomp(topVars)
summary(pca_dat)
# Plot proportions of variances
summary(pca_dat)$importance
summary(pca_dat)$importance[2,]

tiff("screeplot_all_samples.tiff")
barplot(summary(pca_dat)$importance[2,], 
        main="Proportion of variance attributed to PCs",
        ylim=c(0, 1))
dev.off()

# Plot rotations (samples)
rotation <- pca_dat$rotation
rotation <- data.frame(rotation, Group=as.character(colData(dds)$Group),
                       Organ=as.character(colData(dds)$Organ),
                       samples=colData(dds)$SampleID)
head(rotation)

ggplot(rotation, aes(PC1, PC2, color=Group, shape=Organ)) + 
  geom_point(size=4, stroke=2)
ggsave("all_samples_PC1_vs_PC2_top500.tiff", device = "tiff",
        units = "in", width = 6, height = 6)

###############################################################################
## ------------------------ FUNCTIONS ---------------------------------------##

# Function to compare 2 conditions
runPairWise <- function(dds, selected, contrast, Organ, Treated, 
                        Untreated, annot) {
  # dds  - DESeq data set object
  # selected - index of columns to select
  # contrast - character vector: condition, treated, untreated
  # Organ - tissue, organ
  # Treated - treated group
  # Untreated - untreated control
  # annot - annotation data frame
  dds_diff <- dds[,selected]
  dds_diff$Group <- droplevels(dds_diff$Group)
  design(dds_diff) <- formula(~ Group)
  # Remove genes with less than 10 reads in at 
  # least 2 samples
  keep <- which(rowSums(counts(dds_diff) >= 10) > 2) 
  dds_diff <- dds_diff[keep,]
  dds_diff <- DESeq(dds_diff)
  save(dds_diff, file=paste(Organ, Treated, 
                            "vs", Untreated, 
                            "diff.RData", sep="_"))
  res <- results(dds_diff, 
                 contrast = contrast)
  
  # Build MA plot
  tiff(file=paste(Organ, Treated, 
                  "vs", Untreated, 
                  "MAplot.tiff", sep="_"), width=500, height=500)
  DESeq2::plotMA(res, ylim=c(-3, 3), cex = 1)
  dev.off()
  
  # Build a heatmap based on significantly changed genes
  vsd_diff <- vst(dds_diff)
  vsd_sig <- assay(vsd_diff[which(res$padj < 0.05),])
  df <- data.frame(row.names = rownames(colData(vsd_diff)),
                   Group = colData(vsd_diff)[,c("Group")])
  tiff(paste(Organ, Treated, 
             "vs", Untreated, 
             "heatmap.tiff", sep="_"), width=600, 
       height = 600)
  pheatmap(vsd_sig, scale = "row", fontsize_row=3,
           color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "RdYlBu")))(100),
           annotation_col=df, 
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           labels_col = dds_diff$SampleID,
           show_rownames = F,
           main=paste(Organ, ":", Treated, 
                      "vs", Untreated, sep=" "))
  dev.off()
  
  # Build volcano plot
  res <- res[with(res, order(log2FoldChange)),]
  res$threshold <- as.factor(res$padj < 0.05)
  res <- res[!is.na(res$padj),]
  ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj), 
                       colour=threshold)) + 
    geom_point(alpha=0.4, size=1.75) + xlim(c(min(res$log2FoldChange), 
                                              c(max(res$log2FoldChange)))) +
    ylim(c(min(-log10(res$padj)), max(-log10(res$padj)))) + 
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(axis.text=element_text(size=12, face="bold")) +
    theme(axis.title=element_text(size=14)) +
    theme(legend.title=element_text(size=14)) +
    theme(legend.text=element_text(size=12))
  ggsave(paste(Organ, Treated, 
               "vs", Untreated, 
               "volcano_plot.tiff", sep="_"), dpi=300, units = "in",
         device = "tiff", height = 5, width = 5)
  
  # Get normalized counts
  norm_counts <- counts(dds_diff, normalized = T)
  
  # Attach normalized counts to the results data frame
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res_counts <- merge(res, norm_counts, by=0)
  res_counts <- res_counts[order(res_counts$padj, decreasing = F),]
  res_annot <- merge(res_counts, annot, by.x="Row.names",
                     by.y="ensembl_gene_id")
  
  write.csv(res_annot, 
            file=paste(Organ, Treated, 
                       "vs", Untreated, 
                       "results.csv", sep="_"),
            row.names = F)
  
  write.csv(res_annot[res_annot$padj < 0.05,], 
            file=paste(Organ, Treated, 
                       "vs", Untreated, 
                       "sig_only.csv", sep="_"),
            row.names = F)
}

############################################################################
## --------------- RUN PAIR-WISE COMPARISONS -----------------------------##
## Get index of samples of interest using any convenient method
selected <- grep("Some_treatment", subset(colData(dds), Condition == "Some_treatment")$Group)
runPairWise(dds=dds, selected=selected, 
            contrast=c("Condition", 
                                "Treated",
                                "Untreated"), 
            Organ="Brain", Treated="Treated", 
            Untreated="Untreated", annot=annot)


