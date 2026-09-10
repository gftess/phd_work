# =========================================================================== #
#### ===================== PDX VS PATIENT - CITOKYNES ==================== ####
# =========================================================================== #
# Date:  03.set.2026 (Happy Brazilian Biologist's Day!)
### How citokynes expression changes with the transplantation of leukemia
## cells in mice?
# <-- INPUT:
## + normalized BRB-seq expression matrix 
## + samples data 
## + target genes list
# --> OUTPUT:
## + log2 Fold Change matrix (.csv)
# =========================================================================== #
#### 0. Libraries -------------------------------------------------------- ####
# =========================================================================== #
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
# =========================================================================== #
#### 1. Inputs ----------------------------------------------------------- ####
# =========================================================================== #
## 1.1. Normalized Brb-seq expression matrix ------------------------------- ##
lcpm <- readRDS("./brbseq/input/20260902_brbseq_normCount_wihtCombat.RDS")
## 1.2. Samples data ------------------------------------------------------- ##
smp.dt <- data.frame("sample" = colnames(lcpm),
                     "condition" = NA)
smp.dt$id <- sub("_.*", "", smp.dt$sample)
keep.pair <- unique(smp.dt$id[duplicated(smp.dt$id)])
smp.dt <- smp.dt[smp.dt$id %in% keep.pair,]
smp.dt <- smp.dt[!grepl("PDX4", smp.dt$sample),]
smp.dt <- smp.dt[-30,] # fixing 200271 problem
smp.dt <- smp.dt[!grepl("PDXpur", smp.dt$sample),]
rownames(smp.dt) <- NULL
smp.dt[276,1] <- "200271_D0_PDX" # 200271 problem fixed
smp.dt <- smp.dt[!grepl("PDX1", smp.dt$sample),]
smp.dt <- smp.dt[!grepl("patient1", smp.dt$sample),]
smp.dt$id <- as.numeric(smp.dt$id)
smp.dt <- smp.dt[order(smp.dt$id), ]
rownames(smp.dt) <-NULL
temp <- as.data.frame(table(smp.dt$id))
smp.dt <- smp.dt[-c(7,8,21,94:98,54,55,66,67,195,196,237,244),]
## 1.3. Citokyne genes ----------------------------------------------------- ##
ck.gene <- read.csv("brbseq/input/citokyne_genes.csv")
## 1.x. Theme for ggplot2 -------------------------------------------------- ##
theme_set(theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                          vjust = 1.5, size = 14),
                axis.title = element_text(hjust = 0.5, face = "bold",
                                          size = 11),
                axis.text = element_text(color = "black", size = 9),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_rect(color = "black"),
                legend.title = element_text(size = 10, face = "bold")))
## 1.y. Output path -------------------------------------------------------- ##
out <- "./brbseq/output/pdx_patient/"
# =========================================================================== #
#### 2. Calculate log2FoldChange ----------------------------------------- ####
# =========================================================================== #
## 2.1. Create a matrix to each condition ---------------------------------- ##
pct <- smp.dt[grepl("patient", smp.dt$sample),1]
pdx <- smp.dt[grepl("PDX", smp.dt$sample),1]
pct.mtx <- lcpm[, pct]
pdx.mtx <- lcpm[, pdx]
## 2.2. Direct subtraction to generate delta matrix ----------------------- ##
lfc.mtx <- pdx.mtx - pct.mtx
# =========================================================================== #
#### 3. Filter citokyne genes -------------------------------------------- ####
# =========================================================================== #
## 3.1. Convert Symbol to Ensembl ------------------------------------------ ##
ck.gene$ensembl <- mapIds(org.Hs.eg.db, 
                      keys = ck.gene$symbol, 
                      keytype = "SYMBOL", 
                      column = "ENSEMBL")
ck.gene <- unique(ck.gene)
## 3.2. Filter genes in matrix -------------------------------------------- ##
lfc.mtx <- as.data.frame(lfc.mtx)
lfc.mtx <- lfc.mtx[intersect(ck.gene$ensembl, rownames(lfc.mtx)),]
## 3.3. Convert Ensembl in Symbol in log2FoldChange matrix ---------------- ##
lfc.mtx$ensembl <- rownames(lfc.mtx)
temp <- merge(lfc.mtx, ck.gene, by = "ensembl", all.x = T)
rownames(temp) <- temp$symbol
temp$ensembl <- NULL
temp$symbol <- NULL
## 3.4. Save log2FoldChange matrix in .csv format ------------------------ ##
write.csv(temp, "./brbseq/output/pdx_patient/20260903_log2FoldChange_matrix.csv")





