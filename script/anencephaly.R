# =========================================================================== #
#### ============== ANENCEPHALY - MICROARRAY RNA EXPRESSION ============== ####
# =========================================================================== #
# Date: 24.ago.2026
### Using expression data from a microRNA microarray assay to perform DEG
## analysis between anencephaly blood samples and healthy controls. 
# <-- INPUT:
## + normalized log2 intensities matrix
## + samples data
# --> OUTPUT:
## + distribution plots (.png)
## + histogram and elbow plot of sd (.png)
## + pca biplots (.png)
## + 
# =========================================================================== #
#### 0. Libraries -------------------------------------------------------- ####
# =========================================================================== #
library(ggplot2)
library(readr)
library(limma)
library(tidyverse)
library(psych)
library(patchwork)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
# =========================================================================== #
#### 1. Inputs ----------------------------------------------------------- ####
# =========================================================================== #
## 1.1. Normalized log2 intensities matrix with sample data ---------------- ##
int.mtx <- read_tsv("./anencephaly/input/intensidades_A_hsa.txt")
int.mtx <- as.data.frame(int.mtx)
colnames(int.mtx) <- c("ID", "AN1", "AN2", "AN3", "AN4", "AN5", "AN6", 
                       "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6")
## 1.2. Theme for ggplot2 -------------------------------------------------- ##
theme_set(theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                          vjust = 1.5, size = 14),
                axis.title = element_text(hjust = 0.5, face = "bold",
                                          size = 11),
                axis.text = element_text(color = "black", size = 9),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_rect(color = "black"),
                legend.title = element_text(size = 10, face = "bold")))
## 1.3. Output path -------------------------------------------------------- ##
out <- "./anencephaly/output/"
# =========================================================================== #
#### 2. Design data ------------------------------------------------------ ####
# =========================================================================== #
data <- data.frame("sample" = colnames(int.mtx),
                   "condition" = ifelse(grepl("CTRL", colnames(int.mtx)),
                                        "Control","Anencephaly"))
rownames(data) <- data$sample
data <- data[-1,]
# =========================================================================== #
#### 3. QC analysis ------------------------------------------------------ ####
# =========================================================================== #
# --------------------------------------------------------------------------- #
# |   Because this data has already been processed and normalized,          | #
# | this QC analysis aims to assess the quality of normalization and the    | #
# | distribution of processed intensities.                                  | #
# --------------------------------------------------------------------------- #
## 3.1. Create a adequate object for ggplot visualization ------------------ ##
int.long <- int.mtx %>% pivot_longer(cols = -ID,
                                     names_to = "sample",
                                     values_to = "count")
gg.dt <- merge(int.long, data, by = "sample", all.x = T)
gg.dt$count <- as.numeric(gsub(",", ".", as.character(gg.dt$count))) 

## 3.2. Distribution of intensity by sample -------------------------------- ##
ggplot(gg.dt,aes(x = sample, y = count, fill = condition)) +
  geom_boxplot(color = "black", outliers = T)+
  ylab("Log2 Intensity")+
  xlab("Sample")+
  labs(title = "Distribution of intensity by sample")+
  guides(fill=guide_legend(title="Condition"))+
  scale_fill_manual(values =c("brown2", "darkslategray3"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("260909_03_qc_distribution_intensity_bySample.png",
       path = out, width = 7, height = 5)
dev.off()
## 3.3. Distribution of overall intensity ---------------------------------- ##
ggplot(gg.dt, aes(x = count)) +
  geom_density(alpha = 0.3, fill = "plum4") +
  labs(title = "Density plot of overall intensity",
       x = "log2 Intensity", y = "Density")
ggsave("260909_03_qc_density_overall_intensity.png",
       path = out, width = 7, height = 5)
dev.off()
## 3.4. Distribution of intensity by sample -------------------------------- ##
ggplot(gg.dt, aes(x = count, group = sample, fill = condition)) +
  geom_density(alpha = 0.3) +
  labs(title = "Density plot intensity by sample",
       x = "log2 Intensity", y = "Density")+
  guides(fill=guide_legend(title="Condition"))+
  scale_fill_manual(values =c("brown2", "darkslategray3"))
ggsave("260909_03_qc_density_bySample_intensity.png",
       path = out, width = 7, height = 5)
dev.off()
# =========================================================================== #
#### 4. Principal Component Analysis (PCA) ------------------------------- ####
# =========================================================================== #
## 4.1. In intensity matrix, numeric values only --------------------------- ##
rownames(int.mtx) <- int.mtx$ID
int.mtx$ID <- NULL
rnames <- rownames(int.mtx)
cnames <- colnames(int.mtx)
int.mtx <- as.matrix(int.mtx)
int.mtx <- matrix(as.numeric(gsub(",", ".", int.mtx)),
                  nrow = nrow(int.mtx), 
                  ncol = ncol(int.mtx))
int.mtx <- as.data.frame(int.mtx)
rownames(int.mtx) <- rnames
colnames(int.mtx) <- cnames
## 4.2. Using all transcripts ---------------------------------------------- ##
## 4.2.1. Defining PCs ----------------------------------------------------- ##
pc <- prcomp(t(int.mtx),
             center = TRUE,
             scale. = TRUE)
summary(pc) 
pairs.panels(pc$x[,1:5],
             gap=0, 
             pch=20)
pca.dt <- data.frame(PC1 = pc$x[, 1],
                     PC2 = pc$x[, 2],
                     PC3 = pc$x[, 3],
                     PC4 = pc$x[, 4],
                     PC5 = pc$x[, 5])
## 4.2.2. Merge with samples data ------------------------------------------ ##
pca.dt$sample <- rownames(pca.dt)
pca.m <- merge(pca.dt, data, by = "sample")
## 4.2.3. Plotting PCs ----------------------------------------------------- ##
pcI <- "PC1 (21.4%)"
pcII <- "PC2 (16.8%)"
pcIII <- "PC3 (12.5%)"
pcIV <- "PC4 (10.3%)"
pcV <- "PC5 (7.7%)"

ggI.II <- ggplot(pca.m, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC2",
       x = pcI,
       y = pcII) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))
ggI.III <- ggplot(pca.m, aes(x = PC1, y = PC3, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC3",
       x = pcI,
       y = pcIII) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))
ggI.IV <- ggplot(pca.m, aes(x = PC1, y = PC4, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC4",
       x = pcI,
       y = pcIV) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(plot.title = element_text(size = 12))
ggI.V <- ggplot(pca.m, aes(x = PC1, y = PC5, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC5",
       x = pcI,
       y = pcV) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(plot.title = element_text(size = 12))

(ggI.II + ggI.III + ggI.IV + ggI.V) +
  plot_annotation(title = "PCA for all transcripts")

ggsave("260909_04_pca_panel_all_transcripts.png",
       path = out, width = 7, height = 7)
dev.off()
## 4.3. Using Most Variable Transcripts (MVT) ------------------------------ ##
## 4.3.1. Visualize the overall standard deviation ------------------------- ##
sd <- apply(int.mtx,1,sd)
png(file.path(out, "260909_04_histogram_sd_all_transcripts.png"), 
    width = 6, height = 4, units = "in", res = 300)
hist(sd, breaks = 100, col = "snow",
     xlab = "Standard Deviation", 
     ylab = "Frequency", main = "Histogram of standard deviation")
abline(v=quantile(sd, 0.50), col = "gray10",lwd = 3, lty = 1)#50%
abline(v=quantile(sd, 0.75), col = "gray30",lwd = 3, lty = 1)#75%
abline(v=quantile(sd, 0.90), col = "gray50",lwd = 3, lty = 1)#90%
abline(v=quantile(sd, 0.99) , col = "gray70",lwd = 3, lty = 1)#95%
legend("right", inset=.02, title="Quantile",
       c("50%","25%","10%","1%"),
       fill=c("gray10", "gray30", "gray50", "gray70"),
       horiz=F, cex=0.8)
box(lwd = 1.5)
dev.off()
## 4.3.2. Elbow plot ------------------------------------------------------- ##
# --------------------------------------------------------------------------- #
# | The inflection point ("elbow") of the curve indicates the threshold     | #
# | where the variation ceases to represent a relevant biological signal    | #
# | and becomes merely background noise.                                    | #
# --------------------------------------------------------------------------- #
## A data frame ordered by highest sd to lowest ##
df.sd <- data.frame(miRNA = names(sd),
                    SD = sd) %>%
  arrange(desc(SD)) %>%
  mutate(Rank = row_number())

ggplot(df.sd, aes(x = Rank, y = SD)) +
  geom_line(color = "paleturquoise2", size = 0.75) +
  geom_point(color = "slateblue", alpha = 0.3, size = 1.5) +
  geom_vline(xintercept = c(100, 250, 500, 1000, 2000), linetype = "dashed",
             color = "orangered", alpha = 0.7) +
  labs(
    title = "Elbow plot of miRNAs by standard deviation",
    x = "miRNAs rank (highest to lowest sd)",
    y = "Standard deviation ")
ggsave("260909_04_elbow_plot.png",
       path = out, width = 7, height = 5)
dev.off()
## 4.3.3. Filter MVT ------------------------------------------------------- ## 
mvt.mtx <- int.mtx[names(sort(sd, decreasing = T))[1:500],]
## 4.3.4. Defining PCs ----------------------------------------------------- ##
pc <- prcomp(t(mvt.mtx),
             center = TRUE,
             scale. = TRUE)
summary(pc) 
pairs.panels(pc$x[,1:5],
             gap=0, 
             pch=20)
pca.dt <- data.frame(PC1 = pc$x[, 1],
                     PC2 = pc$x[, 2],
                     PC3 = pc$x[, 3],
                     PC4 = pc$x[, 4],
                     PC5 = pc$x[, 5])
## 4.3. Merge with samples data -------------------------------------------- ##
pca.dt$sample <- rownames(pca.dt)
pca.m <- merge(pca.dt, data, by = "sample")
## 4.4. Plotting PCs ------------------------------------------------------- ##
pcI <- "PC1 (27.5%)"
pcII <- "PC2 (20.9%)"
pcIII <- "PC3 (13.3%)"
pcIV <- "PC4 (9.8%)"
pcV <- "PC5 (6.4%)"

ggI.II <- ggplot(pca.m, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC2",
       x = pcI,
       y = pcII) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))
ggI.III <- ggplot(pca.m, aes(x = PC1, y = PC3, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC3",
       x = pcI,
       y = pcIII) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))
ggI.IV <- ggplot(pca.m, aes(x = PC1, y = PC4, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC4",
       x = pcI,
       y = pcIV) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(plot.title = element_text(size = 12))
ggI.V <- ggplot(pca.m, aes(x = PC1, y = PC5, color = condition)) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 2, color = "black")+
  labs(title = "PC1 x PC5",
       x = pcI,
       y = pcV) +
  guides(color = "none")+
  scale_color_manual(values = c("brown2", "darkslategray3"))+
  theme(plot.title = element_text(size = 12))

(ggI.II + ggI.III + ggI.IV + ggI.V) +
  plot_annotation(title = "PCA for the 100 most variable miRNAs")

ggsave("260909_04_pca_panel_100_mvt.png",
       path = out, width = 7, height = 7)
dev.off()
# =========================================================================== #
#### 5. Limma and linear model ------------------------------------------- ####
# =========================================================================== #
## 5.0. Remove anencephaly samples that were clustered with control samples  ##
int.mtx <- int.mtx[,c("AN1","AN6", "AN5", 
                      "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6")]
data <- data[data$sample %in% c("AN1","AN6", "AN5", 
                                "CTRL1", "CTRL2", "CTRL3", "CTRL4",
                                "CTRL5", "CTRL6"),]
## 5.1. Make sure the data is a numeric matrix ----------------------------- ##
int.mtx <- as.matrix(int.mtx)
## 5.2. Create design ------------------------------------------------------ ##
data <- data[colnames(int.mtx),]
condition <- factor(data$condition)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
## 5.3. Create contrast matrix --------------------------------------------- ##
contrast.mtx <- makeContrasts(Anencephaly_vs_Control = Anencephaly - Control,
                           levels = design)
## 5.4. Adjust linear model ------------------------------------------------ ##
fit <- lmFit(int.mtx, design)
## 5.5. Apply contrasts ---------------------------------------------------- ##
fit2 <- contrasts.fit(fit, contrast.mtx)
## 5.6. Moderate standard errors with Empirical Bayes (eBayes) ------------- ##
fit2 <- eBayes(fit2)
## 5.7. Results ------------------------------------------------------------ ##
res <- topTable(fit2, coef = "Anencephaly_vs_Control",
                              number = Inf, adjust.method = "BH")
deg.res <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.58) ## 0.58 = alteration of 1.5x
## 5.4. Define miRNA set to be analyzed (all | high | low in AN) ----------- ##
all <- deg.res
low <- deg.res[deg.res$logFC < 0,]
high <- deg.res[deg.res$logFC > 0,]
# =========================================================================== #
#### 6. Enrichment Analysis ---------------------------------------------- ####
# =========================================================================== #
# --------------------------------------------------------------------------- #
# | Make sure the miRNA names are the official ones, such as those with     | #
# | 'hsa' and similar. The transformation needed for Affymetrix array are   | #
# | performed in this script.                                               | #
# --------------------------------------------------------------------------- #
## 6.1. Convert names to multiMiR compatible format ------------------------ ##
my.mirnas <- rownames(high) # <<<< DEFINE SET DE miRNAs HERE <<<<
my.mirnas <- sub("_st$", "", my.mirnas)
my.mirnas <- my.mirnas[grep("^hsa-", my.mirnas)] ## << human
my.mirnas <- unique(my.mirnas)
## 6.2. Map validated target genes on multiMiR ----------------------------- ##
res.multimir <- get_multimir(mirna = my.mirnas,
                             # consult miRTarBase, TarBase e BioGRID
                             table = "validated",
                             summary = TRUE)
targets.df <- res.multimir@data
## 6.3. Select targets with strong support --------------------------------- ##
targets.strg <- targets.df %>%
  filter(support_type == "Functional MTI" | is.na(support_type))
## 6.4. Extract unique Entrez ID ------------------------------------------- ##
targets.entrez <- unique(targets.strg$target_entrez)
targets.entrez <- targets.entrez[!is.na(targets.entrez)]
## 6.5. Enrichment of Gene Ontology (GO) terms ----------------------------- ##
ego <- enrichGO(gene = targets.entrez,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",  # Biological Process
                pAdjustMethod = "BH",  # Benjamini-Hochberg correction = FDR
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
## 6.5.1. Reduce redundancy of very similar GO terms ----------------------- ##
ego.s <- simplify(ego, 
                  cutoff = 0.7, 
                  by = "p.adjust", 
                  select_fun = min)
ego.s <- setReadable(ego.s,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID")
## 6.6. Enrichment of Kegg pathways ---------------------------------------- ##
ekegg <- enrichKEGG(gene = targets.entrez,
                    organism = "hsa",
                    pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID")
## 6.7. Visualize ---------------------------------------------------------- ##
## 6.7.1. GO --------------------------------------------------------------- ##
## 6.7.1.1. Classical Dotplot ---------------------------------------------- ##
p1 <- dotplot(ego.s, showCategory = 10) + 
  ggtitle("GO termns (Biological Process)")
print(p1)
## 6.7.1.2. Cnetplot ------------------------------------------------------- ##
p2 <- cnetplot(ego.s, categorySize = ~pvalue, foldChange = NULL) +
  ggtitle("Network of Biological Process Terms")
print(p2)

(p1+p2)+
  plot_annotation(title = "GO Enrichment Analysis for low DE miRNAs")+
  plot_layout(widths = c(1, 2.5))
ggsave("260909_06_enrich_highDE.png",
       path = out, width = 15, height = 6)
dev.off()
## 6.7.2. KEGG ------------------------------------------------------------- ##
## 6.7.2.1. Classical Dotplot ---------------------------------------------- ##
p3 <- dotplot(ekegg, showCategory = 10) + 
  ggtitle("KEGG Pathways")
print(p3)
## 6.7.2.2. Cnetplot ------------------------------------------------------- ##
p4 <- cnetplot(ekegg, categorySize = ~pvalue, foldChange = NULL) +
  ggtitle("Network of KEGG Pathways")
print(p4)

(p3+p4)+
  plot_annotation(title = "KEGG Enrichment Analysis for low DE miRNAs")+
  plot_layout(widths = c(1, 2.5))
ggsave("260909_06_enrich_KEGG_highDE.png",
       path = out, width = 15, height = 6)
dev.off()

