# =========================================================================== #
#### =========================== PRONON - GSEA =========================== ####
# =========================================================================== #
# Date: 18.ago.2026
### GSEA analysis using auc matrix with results of ELDA drug screening and 
## BRB-seq data of the same samples. 
# --> INPUT:
## + auc matrix
## + samples brb-seq batch
## + brb-seq raw count matrix
## + analysis samples set
# --> OUTPUT:
## + 
# =========================================================================== #
#### 0. Libraries -------------------------------------------------------- ####
# =========================================================================== #
library(ggplot2)
library(DESeq2)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(ggrepel)
# =========================================================================== #
#### 1. Inputs ----------------------------------------------------------- ####
# =========================================================================== #
## 1.1. Drug screening AUC matrix ##
# It has been already filtered by sample set, using the dashboard 
# "Download Tables" 
auc.mtx <- read.csv("./pronon_gsea/input/auc-matrix.csv")
#auc.mtx$plate_1 <- NULL
#auc.mtx$plate_2 <- NULL

## 1.2. Samples brb-seq batch ##
batch.brb <- read.csv("../msc_work/alithea_brbseq/input/20260818_alithea_batch.csv")
batch.brb$X <- NULL

## 1.3. BRB-seq raw count matrix ##
brb.mtx <- readRDS("../msc_work/alithea_brbseq/input/20260824_brbseq_umi_count_raw_MQUDI1to6.RDS")

## 1.4. Analysis samples set ##
smp.set <- read.csv("pronon_gsea/input/analysis_sample_set_BALL_60h.csv")
smp.set <- smp.set[smp.set$brbseq != "sus",] # samples with problems
smp.set <- smp.set[smp.set$brbseq != "",] # without data
ss.plate1 <- na.omit(smp.set$plate_1)
ss.plate2 <- na.omit(smp.set$plate_2)
smp.set$brbseq <- sub(" ", "", smp.set$brbseq)
# ------------------------ #
# | >>>> 85 SAMPLES <<<< | #
# ------------------------ #
# =========================================================================== #
#### 2. Data cleaning ---------------------------------------------------- ####
# =========================================================================== #
## 2.1. Filter samples set in data ##
#colnames(brb.mtx) <- sub("X", "", colnames(brb.mtx))
smp.set$brbseq <- sub(" ", "", smp.set$brbseq)
brb.mtx <- brb.mtx[,smp.set$brbseq]

#auc.mtx <- auc.mtx[auc.mtx$sample %in% smp.set$ext_timepoint_name,]
#temp <- subset(smp.set, select = c(ext_timepoint_name, brbseq))
#temp <- merge(temp, auc.mtx, by.x = "ext_timepoint_name", by.y = "sample")
#rownames(temp) <- temp$brbseq
#temp$ext_timepoint_name <- NULL
#temp$brbseq <- NULL
#auc.mtx <- temp
# =========================================================================== #
#### 3. VST -------------------------------------------------------------- ####
# =========================================================================== #
# --------------------------------------------------------------------------- #
# |   This part is adapted from run_vst.R script, kept in 2025_ELDA_Omics   | #
# | repository in gitlab.unicamp and developed by PhD. Guilherme Giusti and | #
# | Professor João Meidanis.                                                | #
# --------------------------------------------------------------------------- #
## 3.1. Convert to numeric matrix for DESeq2 ##
counts <- as.matrix(brb.mtx)
storage.mode(counts) <- "numeric"

## 3.2. QC plots of raw counts ##
plot_hist <- function(x, png_path, main, xlab, bins = 80, use_log10 = FALSE) {
  png(png_path, width = 1800, height = 1200, res = 300)
  if (use_log10) {
    x <- log10(x + 1)
  }
  hist(x, breaks = bins, main = main, 
       xlab = xlab, ylab = "Frequency", border = NA)
  dev.off()
}

gene_sums   <- rowSums(counts, na.rm = TRUE)  # total UMI per gene
sample_sums <- colSums(counts, na.rm = TRUE)  # total UMI per sample

plot_hist(gene_sums,
          file.path("./", "gene_umi_hist.png"),
          "Per-gene total UMI distribution",
          "Total UMI per gene",
          bins = 80,
          use_log10 = FALSE)

plot_hist(gene_sums,
          file.path("./", "gene_umi_hist_log.png"),
          "Per-gene total UMI distribution (log10)",
          "log10(total UMI per gene + 1)",
          bins = 80,
          use_log10 = TRUE)

plot_hist(sample_sums,
          file.path("./", "sample_umi_hist.png"),
          "Per-sample library size distribution",
          "Total UMI per sample",
          bins = 30,
          use_log10 = FALSE)

plot_hist(sample_sums,
          file.path("./", "sample_umi_hist_log.png"),
          "Per-sample library size distribution (log10)",
          "log10(total UMI per sample + 1)",
          bins = 30,
          use_log10 = TRUE)

## 3.3. Build DESeq2 object with neutral design ##
# I changed to correct the batch effect ##
coldata <- batch.brb[batch.brb$PSID %in% colnames(brb.mtx),]
rownames(coldata) <- NULL
coldata <- coldata[-18,]
rownames(coldata) <- NULL
coldata <- coldata[-4,]
rownames(coldata) <- NULL
coldata <- coldata[-5,]
rownames(coldata) <- coldata$PSID
coldata <- coldata[colnames(brb.mtx),]
coldata$PSID <- NULL

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData   = coldata,
                              design    = ~ Group)

min_count   <- 5
min_samples <- 3
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds  <- dds[keep, ]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

## 3.3. Normalized Counts ##
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid")
write.csv(norm_counts, file.path("./pronon_gsea", "normalizedcounts.csv"),
          row.names = FALSE)

## 3.4. Variance Stabilizing Transformation ##
vsd <- vst(dds, blind = FALSE) # changed to FALSE because of batch effect
vst_mat <- limma::removeBatchEffect(assay(vsd), batch = coldata$Group) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Geneid")
write.csv(vst_mat, file.path("./pronon_gsea", "vst_counts.csv"),
          row.names = FALSE)




