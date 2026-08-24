# =========================================================================== #
#### =========================== ALL SECRETOME =========================== ####
# =========================================================================== #

# Date: 24.abr.2026
### How many secreted proteins of the human genome are expressed in ALL? 

# --> INPUT:
## + list of secreted proteins from uniprot
## + alithea brb-seq normlized count matriz (log2CPM)
## + alithea metadata

# --> OUTPUT:
## + .csv files with filtered genes
## + .png file with Venn Diagram
## + .png files with enrichment analysis

# =========================================================================== #
#### 0. Libraries -------------------------------------------------------- ####
# =========================================================================== #

library(ggplot2)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggVennDiagram)

# =========================================================================== #
#### 1. Input ------------------------------------------------------------ ####
# =========================================================================== #

## 1.1. List of secreted proteins ####
plist <- read.csv("input/uniprot_secreted.csv",
                  sep = ",")
plist <- plist %>%
  separate_rows(Gene.Names, sep = " ")

## 1.2. BRB-seq ####
## 1.2.1. Normalized count matrix ##
#brb.c <- readRDS("../../msc_work/alithea_brbseq/output/
#20260424_alithea_count_RLE_log2CPM_wMQUDI5.rds")
## 1.2.2. Metadata alithea for subtype information ##
brb.ss <- 
  read.csv("../../msc_work/alithea_brbseq/input/alithea_metadata_1_5.csv")
brb.ss <- brb.ss[brb.ss$PSID %in% colnames(brb.c),]
brb.ss <- brb.ss[!duplicated(brb.ss$PSID),]

brb.c <-  brb.c[,c(brb.ss$PSID, "symbol")] ## >>>> 279 amostras ##

# =========================================================================== #
#### 2. Analysis --------------------------------------------------------- ####
# =========================================================================== #

## 2.1. Filter secreted proteins in count matrix ####
brb.c <- brb.c[brb.c$symbol %in% plist$Gene.Names,] ## >>>> 2896 genes ##

## 2.2. Filter leukemia by lineage ####
## 2.2.1. B-ALL ####
ssb <- brb.ss[brb.ss$Subtype != "LLA-T",]
ssb <- ssb[ssb$Subtype != "LLA-T/ABL1-break",]
ssb <- ssb[ssb$Subtype != "LLA T",]

## 2.2.1.1. Patient samples ##
ssbp <- ssb[grep("pcte",ssb$Sample_ID_fixed),]
cbp <- brb.c[, c("symbol", ssbp$PSID)]
symbol <- cbp$symbol
cbp$symbol <- NULL

keep <- as.character()
out <- as.character()
for (a in rownames(cbp)) {
  soma <- sum(cbp[a,]>1)
  if (soma > 10){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
cbp <- cbp[keep,]

keep <- as.character()
out <- as.character()
for (a in rownames(cbp)) {
  z <- t(cbp[a,])
  median <- median(z[,1])
  if (median > 2){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
cbp <- cbp[keep,]

cbp$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(cbp),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

write.csv(cbp, "output/20260427_secretome_BALL_pct.csv")

## 2.2.1.2. PDX samples ##
ssbx <- ssb[grep("PDX",ssb$Sample_ID_fixed),]
cbx <- brb.c[, c("symbol", ssbx$PSID)]

cbx$symbol <- NULL

keep <- as.character()
out <- as.character()
for (a in rownames(cbx)) {
  soma <- sum(cbx[a,]>1)
  if (soma > 10){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
cbx <- cbx[keep,]

keep <- as.character()
out <- as.character()
for (a in rownames(cbx)) {
  z <- t(cbx[a,])
  median <- median(z[,1])
  if (median > 2){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
cbx <- cbx[keep,]

cbx$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(cbx),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

write.csv(cbx, "output/20260427_secretome_BALL_pdx.csv")

## 2.2.2. T-ALL ####
sst <- brb.ss[brb.ss$Subtype %in% c("LLA-T", "LLA-T/ABL1-break",
                                    "LLA T"),]

## 2.2.1.2. Patient samples ##
sstp <- sst[grep("pcte",sst$Sample_ID_fixed),]
ctp <- brb.c[, c("symbol", sstp$PSID)]
symbol <- ctp$symbol
ctp$symbol <- NULL

keep <- as.character()
out <- as.character()
for (a in rownames(ctp)) {
  soma <- sum(ctp[a,]>1)
  if (soma > 10){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
ctp <- ctp[keep,]

keep <- as.character()
out <- as.character()
for (a in rownames(ctp)) {
  z <- t(ctp[a,])
  median <- median(z[,1])
  if (median > 2){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
ctp <- ctp[keep,]

ctp$symbol <- mapIds(org.Hs.eg.db,
                  keys = rownames(ctp),
                  keytype = "ENSEMBL",
                  column = "SYMBOL")

write.csv(ctp, "output/20260427_secretome_TALL_pct.csv")

## 2.2.2.2. PDX samples ##
sstx <- sst[grep("PDX",sst$Sample_ID_fixed),]
ctx <- brb.c[, c("symbol", sstx$PSID)]

ctx$symbol <- NULL

keep <- as.character()
out <- as.character()
for (a in rownames(ctx)) {
  soma <- sum(ctx[a,]>1)
  if (soma > 10){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
ctx <- ctx[keep,]

keep <- as.character()
out <- as.character()
for (a in rownames(ctx)) {
  z <- t(ctx[a,])
  median <- median(z[,1])
  if (median > 2){
    keep <- c(keep, a)
  } else {
    out <- c(out,a)
  }
}
ctx <- ctx[keep,]

ctx$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(ctx),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

write.csv(ctx, "output/20260427_secretome_TALL_pdx.csv")

## 2.3. Venn Diagram ####

ggVennDiagram(list(rownames(cbp),
                   rownames(cbx),
                   rownames(ctp),
                   rownames(ctx)),
              label_alpha = 0,
              category.names = c("BALL_pct","BALL_pdx",
                                 "TALL_pct","TALL_pdx"),
              set_size = 3.5,
              label = c("count"))+
  scale_fill_continuous(palette = c("snow","firebrick4"))+
  theme(legend.position = "none")
res <- erpa@result
ggplot(head(res, 10), aes(x = -log10(qvalue),
                          y = reorder(Description, -qvalue))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  labs(y = "Reactome Pathway",
       x = "-log10(qvalue)", size = "Count") +
  theme_bw()+
  theme(axis.text = element_text(color = "black"))

# =========================================================================== #
#### 3. Enrichment analysis ---------------------------------------------- ####
# =========================================================================== #

## 3.1. Define 'background'  of secreted proteins ####
spbg <- mapIds(org.Hs.eg.db,
                     keys = plist$Gene.Names,
                     keytype = "SYMBOL",
                     column = "ENSEMBL")

## 3.2. cluster profiler - GO (MF, BP, CC) ####
ego <- enrichGO(gene          = rownames(ctp),
                OrgDb         = 'org.Hs.eg.db',
                keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                ont = "CC",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1)

barplot(ego, showCategory=10)

## 3.3. ReactomePA ####
geneid <- mapIds(org.Hs.eg.db,
                 keys = rownames(ctp),
                 keytype = "ENSEMBL",
                 column = "ENTREZID")

erpa <- enrichPathway(geneid,
                      organism = "human",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.1,
                      minGSSize = 10,
                      maxGSSize = 500,
                      readable = FALSE)

# *************************************************************************** #