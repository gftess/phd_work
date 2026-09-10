# =========================================================================== #
#### =========================== PRONON - GSEA =========================== ####
# =========================================================================== #
# Date: 18.ago.2026
### GSEA analysis using auc matrix with results of ELDA drug screening and 
## BRB-seq data of the same samples. 
# --> INPUT:
## + auc matrix or ld50 matrix 
## + samples brb-seq batch
## + brb-seq raw count matrix
## + analysis samples set
## + manual group settings
## + gene map (ENSEMBL to SYMBOL)
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
## 1.1. Drug screening AUC matrix ------------------------------------------ ##
# It has been already filtered by sample set, using the dashboard 
# "Download Tables" 
auc.mtx <- read.csv("./pronon_gsea/auc-matrix.csv")
auc.mtx$plate_1 <- NULL
auc.mtx$plate_2 <- NULL
## 1.2. Samples brb-seq batch ---------------------------------------------- ##
batch.brb <- read.csv("../msc_work/alithea_brbseq/input/20260818_alithea_batch.csv")
batch.brb$X <- NULL
batch.brb[370,1] <- "200271_D0_PDX1"
## 1.3. BRB-seq raw count matrix ------------------------------------------- ##
brb.mtx <- readRDS("../msc_work/alithea_brbseq/input/20260824_brbseq_umi_count_raw_MQUDI1to6.RDS")
## 1.4. Analysis samples set ----------------------------------------------- ##
smp.set <- read.csv("pronon_gsea/analysis_sample_set_TALL_60h.csv")
smp.set <- smp.set[smp.set$brbseq != "sus",] # samples with problems
smp.set <- smp.set[smp.set$brbseq != "",] # without data
ss.plate1 <- na.omit(smp.set$plate_1)
ss.plate2 <- na.omit(smp.set$plate_2)
smp.set$brbseq <- sub(" ", "", smp.set$brbseq)
## 1.5. Manual group settings ---------------------------------------------- ##
man.grp <- read.csv("./pronon_gsea/manual_groups.csv")
## 1.6. Gene map ----------------------------------------------------------- ##
gn.map <- read.csv("./pronon_gsea/gene_map.csv")
gn.map$X <- NULL
# ------------------------ #
# | >>>> 58 SAMPLES <<<< | #
# ------------------------ #
# =========================================================================== #
#### 2. Data cleaning ---------------------------------------------------- ####
# =========================================================================== #
## 2.1. Filter samples set in data ----------------------------------------- ##
#colnames(brb.mtx) <- sub("X", "", colnames(brb.mtx))
smp.set$brbseq <- sub(" ", "", smp.set$brbseq)
brb.mtx <- brb.mtx[,smp.set$brbseq]
#brb.mtx <- brb.mtx[,man.grp$sample]

auc.mtx <- auc.mtx[auc.mtx$sample %in% smp.set$ext_timepoint_name,]
temp <- subset(smp.set, select = c(ext_timepoint_name, brbseq))
temp <- merge(temp, auc.mtx, by.x = "ext_timepoint_name", by.y = "sample")
rownames(temp) <- temp$brbseq
temp$ext_timepoint_name <- NULL
temp$brbseq <- NULL
auc.mtx <- temp
#auc.mtx <- auc.mtx[auc.mtx$sample %in% smp.set$sample,]
auc.mtx$sample <- rownames(auc.mtx)
write.csv(auc.mtx, "./pronon_gsea/auc-matrix.csv", row.names = F)
     # =========================================================================== #
#### 3. VST -------------------------------------------------------------- ####
# =========================================================================== #
# --------------------------------------------------------------------------- #
# |   This part is adapted from run_vst.R script, kept in 2025_ELDA_Omics   | #
# | repository in gitlab.unicamp and developed by PhD. Guilherme Giusti and | #
# | Professor João Meidanis.                                                | #
# --------------------------------------------------------------------------- #
## 3.1. Convert to numeric matrix for DESeq2 --------------------------------##
counts <- as.matrix(brb.mtx)
storage.mode(counts) <- "numeric"
## 3.2. QC plots of raw counts --------------------------------------------- ##
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
          file.path("./pronon_gsea", "gene_umi_hist.png"),
          "Per-gene total UMI distribution",
          "Total UMI per gene",
          bins = 80,
          use_log10 = FALSE)

plot_hist(gene_sums,
          file.path("./pronon_gsea", "gene_umi_hist_log.png"),
          "Per-gene total UMI distribution (log10)",
          "log10(total UMI per gene + 1)",
          bins = 80,
          use_log10 = TRUE)

plot_hist(sample_sums,
          file.path("./pronon_gsea", "sample_umi_hist.png"),
          "Per-sample library size distribution",
          "Total UMI per sample",
          bins = 30,
          use_log10 = FALSE)

plot_hist(sample_sums,
          file.path("./pronon_gsea", "sample_umi_hist_log.png"),
          "Per-sample library size distribution (log10)",
          "log10(total UMI per sample + 1)",
          bins = 30,
          use_log10 = TRUE)
## 3.3. Build DESeq2 object with neutral design ---------------------------- ##
# I changed to correct the batch effect ##
coldata <- batch.brb[batch.brb$PSID %in% colnames(brb.mtx),]
rownames(coldata) <- NULL
coldata <- coldata[-1,]
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
## 3.3. Normalized Counts -------------------------------------------------- ##
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid")
write.csv(norm_counts, file.path("./pronon_gsea", "normalizedcounts.csv"),
          row.names = FALSE)
## 3.4. Variance Stabilizing Transformation -------------------------------- ##
vsd <- vst(dds, blind = FALSE) # changed to FALSE because of batch effect
vst_mat <- limma::removeBatchEffect(assay(vsd), batch = coldata$Group) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Geneid")
write.csv(vst_mat, file.path("./pronon_gsea", "vst_counts.csv"),
          row.names = FALSE)
# =========================================================================== #
#### 4. DESeq ------------------------------------------------------------ ####
# =========================================================================== #
# --------------------------------------------------------------------------- #
# |   This part is adapted from run_deseq.R script, kept in 2025_ELDA_Omics | #
# | repository in gitlab.unicamp and developed by PhD. Guilherme Giusti and | #
# | Professor João Meidanis. In this case, works for manual groups only.    | #
# --------------------------------------------------------------------------- #
## 4.1. Build DESeq2 dataset ------------------------------------------------ #
rownames(man.grp) <- man.grp$sample
man.grp <- man.grp[colnames(counts),]
man.grp$sample <- NULL
man.grp$group <- as.factor(man.grp$group)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = man.grp,
                              design = ~ group)

min_count   <- 5
min_samples <- 3
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds  <- dds[keep, ]
dds <- DESeq(dds)
## 4.2. Save DE results ---------------------------------------------------- ##
res <- results(dds, contrast = c("group", "1", "0"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("Geneid") %>%
  left_join(gn.map, by = "Geneid") %>%
  relocate(Symbol, .after = Geneid)
write.csv(res_df, file.path("./pronon_gsea/", "DESeq2Results.csv"),
          row.names = FALSE)
## 4.3. Save significant genes --------------------------------------------- ##
sigs <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) # >> ZERO!
write.csv(sigs, file.path("./pronon_gsea/", "DEGs.csv"), row.names = FALSE)
## 4.4. Reuse VST outputs: VST for PCA, Normalized for heatmap ------------- ##
vst_cohort  <- read.csv("./pronon_gsea/vst_counts.csv", 
                        row.names = 1, check.names = FALSE)
  norm_cohort <- read.csv("./pronon_gsea/normalizedcounts.csv",
                          row.names = 1, check.names = FALSE)
vst_subset  <- vst_cohort[,  colnames(dds), drop = FALSE]
norm_subset <- norm_cohort[, colnames(dds), drop = FALSE]
## 4.5. PCA plot ----------------------------------------------------------- ##
pca <- prcomp(t(vst_subset), scale. = FALSE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))[1:2]
pcaData <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                      group = man.grp[colnames(vst_subset),
                                      "group", drop = TRUE])
ggplot(pcaData, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=3) +
  labs(x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) +
  theme_classic()
ggsave(file.path("./pronon_gsea/" ,"PCAplot.png"),
       width = 7, height = 5, dpi = 120)
## 4.6. MA Plot with labels (top N DEGs by FDR then |LFC| ------------------ ##
plot_df <- res_df %>%
  mutate(is_sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) %>%
  mutate(Symbol = ifelse(is.na(Symbol) | Symbol == "", Geneid, Symbol))
top_n <- 30
label_df <- plot_df %>%
  filter(is_sig) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n)
ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(data = subset(plot_df, is_sig), alpha = 0.8, size = 1.2) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(
    data = label_df, aes(label = Symbol),
    max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0
  ) +
  labs(x = "log10(baseMean + 1)", y = "log2 fold change",
       title = paste0("MA plot (labeled top ", top_n, " DEGs)")) +
  theme_classic()
ggsave(file.path("./pronon_gsea/", "MAplot_labeled.png"),
       width = 8, height = 6, dpi = 120)
## 4.7. Volcano plot (labeled top N DEGs) ---------------------------------- ##
plot_df_v <- plot_df %>%
  mutate(neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)))
top_n <- 30
label_df <- plot_df_v %>%
  filter(is_sig) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n)
ggplot(plot_df_v, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(data = subset(plot_df_v, is_sig), alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Symbol),
    max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0
  ) +
  labs(x = "log2 fold change",
       y = expression(-log[10](FDR)),
       title = paste0("Volcano plot (labeled top ", top_n, " DEGs)")) +
  theme_classic()
ggsave(file.path("./pronon_gsea/", "Volcano_labeled.png"),
       p_volcano_lab, width = 8, height = 6, dpi = 120)
## 4.8. Dispersion Plot ---------------------------------------------------- ##
png(file.path("./pronon_gsea/", "Dispersionplot.png"),
    width = 800, height = 600)
plotDispEsts(dds)
dev.off()
## 4.9. Heatmap ------------------------------------------------------------ ##
top_genes <- sigs$Geneid
if (length(top_genes) >= 2) {
  norm_counts <- norm_subset %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Geneid")
  mat <- norm_counts %>%
    column_to_rownames("Geneid") %>%
    .[top_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  sym_vec <- sigs$Symbol[match(rownames(mat_z), sigs$Geneid)]
  rownames(mat_z) <- ifelse(is.na(sym_vec) | sym_vec == "",
                            rownames(mat_z),
                            paste0(sym_vec, " (", rownames(mat_z), ")"))
  # Column annotation: sample group
  grp <- as.character(man.grp[colnames(mat_z), "group", drop = TRUE])
  ha <- HeatmapAnnotation(
    df = data.frame(Group = factor(grp, levels = c("0","1"))),
    col = list(Group = c("0" = "blue", "1" = "red"))
  ) 
  png(file.path("./pronon_gsea", "heatmap.png"), width = 1000, height = 1200)
  ht <- Heatmap(mat_z, name = "Z-score",
                cluster_columns = TRUE, cluster_rows = TRUE,
                top_annotation = ha)
  draw(ht)
  dev.off()}
