###############################################
###############################################
# INTEGRATION SCRIPT P53 AT PROMOTER  #
###############################################
###############################################

library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(HiCaptuRe)
library(ggpubr)
library(viridis)

# Nutlin1h ----
# Loading data ----
## functional p53 ----
functional_p53 <- read_rds("../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin1h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin1h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)


DE_H3K27ac_Nutlin1h <-  read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")

# Integrate PCHiC and functional p53 at enhancer
nut1h_ints <- load_interactions(file = paste0("../data/ComplementaryData/PCHiC/Nutlin1h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut1h_ints <- annotate_interactions(nut1h_ints,annotation = "../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

source("../scripts/integration_function.R")
link_p <- p53_genes_interactions(ints=nut1h_ints,p53.chip =functional_p53_p, promoter.annot = promoters,at.promoter=T)
link_p <- link_p[!grepl("non-annotated|\\.",link_p$target_gene),]
link_p <- link_p[link_p$target_gene %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
x <- data.frame(target_gene=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_p <- left_join(link_p,x,multiple = "all")
link_p <- link_p[!is.na(link_p$log2FC),]
link_p <- link_p %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_p$log2FC <- NULL

# promoter target gene 
aux <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
aux <- data.frame(promoter_target_gene=as.character(aux),target_gene=aux$ensembl_gene_id)
link_p <- left_join(link_p,aux)

# log2FC H3K27ac target gene
prom <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
merge <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
merge <- unique(data.frame(target_gene=as.character(merge$ensembl_gene_id),log2FC_H3K27ac_promoter_1h_vs_0h=merge$log2FoldChange))
link_p <- left_join(link_p,merge)

# log2FC expr target gene
colnames(DE_expr_Nutlin1h)[7] <- "target_gene"
link_p <- left_join(link_p,DE_expr_Nutlin1h[,c(2,7)])
colnames(link_p)[13] <- "log2FC_expr_1h_vs_0h_ensg"

colnames(DE_expr_Nutlin10h)[7] <- "target_gene"
link_p <- left_join(link_p,DE_expr_Nutlin10h[,c(2,7)])
colnames(link_p)[14] <- "log2FC_expr_10h_vs_0h_ensg"

aux <- unique(link_p %>% group_by(target_gene) %>% summarize(log2FC_H3K27ac_promoter_1h_vs_0h=mean(log2FC_H3K27ac_promoter_1h_vs_0h,rm.na=T),
                                                           log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                           log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                           p53=p53))
link_p$log2FC_H3K27ac_promoter_1h_vs_0h <- NULL
link_p$log2FC_expr_1h_vs_0h_ensg <- NULL
link_p$log2FC_expr_10h_vs_0h_ensg <- NULL
to.write <- unique(left_join(link_p,aux))

write_rds(to.write,"../data/ComplementaryData/PCHiC_integration/proximal_integration_Nutlin1h.Rds",compress = "xz")


# Nutlin10h ----
# Loading data ----
## functional p53 ----
functional_p53 <- read_rds("../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin10h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin10h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)


DE_H3K27ac_Nutlin1h <-  read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")

# Integrate PCHiC and functional p53 at enhancer
nut10h_ints <- load_interactions(file = paste0("../data/ComplementaryData/PCHiC/Nutlin10h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut10h_ints <- annotate_interactions(nut10h_ints,annotation = "../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

source("../scripts/integration_function.R")
link_p <- p53_genes_interactions(ints=nut10h_ints,p53.chip =functional_p53_p, promoter.annot = promoters,at.promoter=T)
link_p <- link_p[!grepl("non-annotated|\\.",link_p$target_gene),]
link_p <- link_p[link_p$target_gene %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
x <- data.frame(target_gene=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_p <- left_join(link_p,x,multiple = "all")
link_p <- link_p[!is.na(link_p$log2FC),]
link_p <- link_p %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_p$log2FC <- NULL

# promoter target gene 
aux <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
aux <- data.frame(promoter_target_gene=as.character(aux),target_gene=aux$ensembl_gene_id)
link_p <- left_join(link_p,aux)

# log2FC H3K27ac target gene
prom <- promoters[promoters$ensembl_gene_id %in% link_p$target_gene,]
merge <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
merge <- unique(data.frame(target_gene=as.character(merge$ensembl_gene_id),log2FC_H3K27ac_promoter_10h_vs_0h=merge$log2FoldChange))
link_p <- left_join(link_p,merge)

# log2FC expr target gene
colnames(DE_expr_Nutlin1h)[7] <- "target_gene"
link_p <- left_join(link_p,DE_expr_Nutlin1h[,c(2,7)])
colnames(link_p)[13] <- "log2FC_expr_1h_vs_0h_ensg"

colnames(DE_expr_Nutlin10h)[7] <- "target_gene"
link_p <- left_join(link_p,DE_expr_Nutlin10h[,c(2,7)])
colnames(link_p)[14] <- "log2FC_expr_10h_vs_0h_ensg"

aux <- unique(link_p %>% group_by(target_gene) %>% summarize(log2FC_H3K27ac_promoter_10h_vs_0h=mean(log2FC_H3K27ac_promoter_10h_vs_0h,rm.na=T),
                                                             log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                             log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                             p53=p53))
link_p$log2FC_H3K27ac_promoter_10h_vs_0h <- NULL
link_p$log2FC_expr_1h_vs_0h_ensg <- NULL
link_p$log2FC_expr_10h_vs_0h_ensg <- NULL
to.write <- unique(left_join(link_p,aux))

write_rds(to.write,"../data/ComplementaryData/PCHiC_integration/proximal_integration_Nutlin10h.Rds",compress = "xz")
