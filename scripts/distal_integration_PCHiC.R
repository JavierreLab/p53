###############################################
###############################################
# INTEGRATION SCRIPT P53 AT ENHANCER  #
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
functional_p53 <- read_rds("../../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin1h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin1h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)

## Functional p53 at enhancer ----
functional_p53_e <- functional_p53[!functional_p53 %in% functional_p53_p,]
enhancers <- read_rds(file = "../../data/ComplementaryData/enhancers.Rds")
enhancers <- enhancers[grepl("Nutlin1h",enhancers$id),]
functional_p53_e <- subsetByOverlaps(functional_p53_e,enhancers)


DE_H3K27ac_Nutlin1h <-  read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")


# Integrate PCHiC and functional p53 at enhancer
nut1h_ints <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin1h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut1h_ints <- annotate_interactions(nut1h_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

source("../../scripts/integration_function.R")
link_e <- p53_genes_interactions(ints=nut1h_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_e <- link_e[!grepl("non-annotated|\\.",link_e$ensg),]
link <- link_e
link <- link[link$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link <- left_join(link,x,multiple = "all")
link <- link[!is.na(link$log2FC),]
link <- link %>% group_by(p53) %>% filter(log2FC==max(log2FC))

# Annotate enhancers for integrated p53
functional_p53_e <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(link$p53),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(enhancers,functional_p53_e)
aux <- data.frame(p53=as.character(merge$functional_p53_e),enhancer=as.character(merge$enhancers))
link <- left_join(link,aux)

# H3K27ac at p53 1h
regions <- makeGRangesFromDataFrame(link %>% separate(p53,c("chr","start","end"),sep = ":|-"))
merge <- mergeByOverlaps(regions,DE_H3K27ac_Nutlin1h)
merge <- data.frame(p53=as.character(merge$regions),log2FC_H3K27ac_p53_1h_vs_0h=merge$log2FoldChange)
link <- left_join(link,merge)


# H3K27ac at promoter 1h
prom <- promoters[promoters$ensembl_gene_id %in% link$ensg,]
merge <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
merge <- data.frame(ensg=merge$ensembl_gene_id,log2FC_H3K27ac_ensg_1h_vs_0h=merge$log2FoldChange)
link <- left_join(link,merge)

# H3K27ac at enhancer 1h
regions <- makeGRangesFromDataFrame(link %>% separate(enhancer,c("chr","start","end"),sep = ":|-"))
merge <- mergeByOverlaps(regions,DE_H3K27ac_Nutlin1h)
merge <- unique(data.frame(enhancer=as.character(merge$regions),log2FC_H3K27ac_enhancer_1h_vs_0h=merge$log2FoldChange))
link <- left_join(link,merge)

link <- left_join(link,DE_expr_Nutlin1h[,c(2,7)])
colnames(link)[13] <- "log2FC_expr_1h_vs_0h_ensg"

link <- left_join(link,DE_expr_Nutlin10h[,c(2,7)])
colnames(link)[14] <- "log2FC_expr_10h_vs_0h_ensg"

gene_names <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")
gene_names <- as.data.frame(gene_names@elementMetadata)
colnames(gene_names)[1] <- "ensg"
link <- left_join(link,gene_names)

aux <- unique(link %>% group_by(p53) %>% summarize(log2FC_H3K27ac_p53_1h_vs_0h=mean(log2FC_H3K27ac_p53_1h_vs_0h,rm.na=T),
                                                   log2FC_H3K27ac_ensg_1h_vs_0h=log2FC_H3K27ac_ensg_1h_vs_0h,
                                                   log2FC_H3K27ac_enhancer_1h_vs_0h=log2FC_H3K27ac_enhancer_1h_vs_0h,
                                                   log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                   log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                   ensg=ensg))
aux <- unique(aux %>% group_by(p53) %>% summarize(log2FC_H3K27ac_p53_1h_vs_0h=log2FC_H3K27ac_p53_1h_vs_0h,
                                                   log2FC_H3K27ac_ensg_1h_vs_0h=log2FC_H3K27ac_ensg_1h_vs_0h,
                                                   log2FC_H3K27ac_enhancer_1h_vs_0h=mean(log2FC_H3K27ac_enhancer_1h_vs_0h,na.rm=T),
                                                   log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                   log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                   ensg=ensg))
aux <- unique(aux %>% group_by(ensg) %>% summarize(log2FC_H3K27ac_p53_1h_vs_0h=log2FC_H3K27ac_p53_1h_vs_0h,
                                                   log2FC_H3K27ac_ensg_1h_vs_0h=mean(log2FC_H3K27ac_ensg_1h_vs_0h,na.rm=T),
                                                   log2FC_H3K27ac_enhancer_1h_vs_0h=log2FC_H3K27ac_enhancer_1h_vs_0h,
                                                   log2FC_expr_1h_vs_0h_ensg=mean(log2FC_expr_1h_vs_0h_ensg,na.rm=T),
                                                   log2FC_expr_10h_vs_0h_ensg=mean(log2FC_expr_10h_vs_0h_ensg,na.rm=T),
                                                   p53=p53,
                                                   ensg=ensg))

link$log2FC<-NULL
link$log2FC_H3K27ac_p53_1h_vs_0h<-NULL
link$log2FC_H3K27ac_ensg_1h_vs_0h<-NULL
link$log2FC_H3K27ac_enhancer_1h_vs_0h<-NULL
link$log2FC_expr_1h_vs_0h_ensg<-NULL
link$log2FC_expr_10h_vs_0h_ensg<-NULL
link$overlap_I <- NULL
link$overlap_II <- NULL
link$int_type <- NULL
link <- unique(link)
to.write <- left_join(link,aux)


write_rds(to.write,"../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1h.Rds",compress = "xz")

# Nutlin10h ----

# Loading data ----
## functional p53 ----
functional_p53 <- read_rds("../../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin10h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin10h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)

## Functional p53 at enhancer ----
functional_p53_e <- functional_p53[!functional_p53 %in% functional_p53_p,]
enhancers <- read_rds(file = "../../data/ComplementaryData/enhancers.Rds")
enhancers <- enhancers[grepl("Nutlin10h",enhancers$id),]
functional_p53_e <- subsetByOverlaps(functional_p53_e,enhancers)


DE_H3K27ac_Nutlin1h <-  read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")

# Integrate PCHiC and functional p53 at enhancer
nut10h_ints <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin10h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut10h_ints <- annotate_interactions(nut10h_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

source("../../scripts/integration_function.R")
link_e <- p53_genes_interactions(ints=nut10h_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_e <- link_e[!grepl("non-annotated|\\.",link_e$ensg),]
link <- link_e
link <- link[link$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link <- left_join(link,x,multiple = "all")
link <- link[!is.na(link$log2FC),]
link <- link %>% group_by(p53) %>% filter(log2FC==max(log2FC))

# Annotate enhancers for integrated p53
functional_p53_e <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(link$p53),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(enhancers,functional_p53_e)
aux <- data.frame(p53=as.character(merge$functional_p53_e),enhancer=as.character(merge$enhancers))
link <- left_join(link,aux)

# H3K27ac at p53 10h
regions <- makeGRangesFromDataFrame(link %>% separate(p53,c("chr","start","end"),sep = ":|-"))
merge <- mergeByOverlaps(regions,DE_H3K27ac_Nutlin10h)
merge <- data.frame(p53=as.character(merge$regions),log2FC_H3K27ac_p53_10h_vs_0h=merge$log2FoldChange)
link <- left_join(link,merge)


# H3K27ac at promoter 10h
prom <- promoters[promoters$ensembl_gene_id %in% link$ensg,]
merge <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
merge <- data.frame(ensg=merge$ensembl_gene_id,log2FC_H3K27ac_ensg_10h_vs_0h=merge$log2FoldChange)
link <- left_join(link,merge)

# H3K27ac at enhancer 10h
regions <- makeGRangesFromDataFrame(link %>% separate(enhancer,c("chr","start","end"),sep = ":|-"))
merge <- mergeByOverlaps(regions,DE_H3K27ac_Nutlin10h)
merge <- unique(data.frame(enhancer=as.character(merge$regions),log2FC_H3K27ac_enhancer_10h_vs_0h=merge$log2FoldChange))
link <- left_join(link,merge)

link <- left_join(link,DE_expr_Nutlin1h[,c(2,7)])
colnames(link)[13] <- "log2FC_expr_1h_vs_0h_ensg"

link <- left_join(link,DE_expr_Nutlin10h[,c(2,7)])
colnames(link)[14] <- "log2FC_expr_10h_vs_0h_ensg"

gene_names <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")
gene_names <- as.data.frame(gene_names@elementMetadata)
colnames(gene_names)[1] <- "ensg"
link <- left_join(link,gene_names)

aux <- unique(link %>% group_by(p53) %>% summarize(log2FC_H3K27ac_p53_10h_vs_0h=mean(log2FC_H3K27ac_p53_10h_vs_0h,rm.na=T),
                                                   log2FC_H3K27ac_ensg_10h_vs_0h=log2FC_H3K27ac_ensg_10h_vs_0h,
                                                   log2FC_H3K27ac_enhancer_10h_vs_0h=log2FC_H3K27ac_enhancer_10h_vs_0h,
                                                   log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                   log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                   ensg=ensg))
aux <- unique(aux %>% group_by(p53) %>% summarize(log2FC_H3K27ac_p53_10h_vs_0h=log2FC_H3K27ac_p53_10h_vs_0h,
                                                  log2FC_H3K27ac_ensg_10h_vs_0h=log2FC_H3K27ac_ensg_10h_vs_0h,
                                                  log2FC_H3K27ac_enhancer_10h_vs_0h=mean(log2FC_H3K27ac_enhancer_10h_vs_0h,na.rm=T),
                                                  log2FC_expr_1h_vs_0h_ensg=log2FC_expr_1h_vs_0h_ensg,
                                                  log2FC_expr_10h_vs_0h_ensg=log2FC_expr_10h_vs_0h_ensg,
                                                  ensg=ensg))
aux <- unique(aux %>% group_by(ensg) %>% summarize(log2FC_H3K27ac_p53_10h_vs_0h=log2FC_H3K27ac_p53_10h_vs_0h,
                                                   log2FC_H3K27ac_ensg_10h_vs_0h=mean(log2FC_H3K27ac_ensg_10h_vs_0h,na.rm=T),
                                                   log2FC_H3K27ac_enhancer_10h_vs_0h=log2FC_H3K27ac_enhancer_10h_vs_0h,
                                                   log2FC_expr_1h_vs_0h_ensg=mean(log2FC_expr_1h_vs_0h_ensg,na.rm=T),
                                                   log2FC_expr_10h_vs_0h_ensg=mean(log2FC_expr_10h_vs_0h_ensg,na.rm=T),
                                                   p53=p53,
                                                   ensg=ensg))

link$log2FC<-NULL
link$log2FC_H3K27ac_p53_1h_vs_0h<-NULL
link$log2FC_H3K27ac_ensg_1h_vs_0h<-NULL
link$log2FC_H3K27ac_enhancer_1h_vs_0h<-NULL
link$log2FC_expr_1h_vs_0h_ensg<-NULL
link$log2FC_expr_10h_vs_0h_ensg<-NULL
link$overlap_I <- NULL
link$overlap_II <- NULL
link$int_type <- NULL
link <- unique(link)
to.write <- left_join(link,aux)


write_rds(to.write,"../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10h.Rds",compress = "xz")




## Distal integration DYNAMISM Nutlin1h vs. DMSO - Figure 6 ----

## functional p53 ----
functional_p53 <- read_rds("../../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin1h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin1h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)

## Functional p53 at enhancer ----
functional_p53_e <- functional_p53[!functional_p53 %in% functional_p53_p,]
enhancers <- read_rds(file = "../../data/ComplementaryData/enhancers.Rds")
enhancers <- enhancers[grepl("Nutlin1h",enhancers$id),]
functional_p53_e <- subsetByOverlaps(functional_p53_e,enhancers)


DE_H3K27ac_Nutlin1h <-  read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")


## PCHi-C Dynamism ----
dmso <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
dmso <- annotate_interactions(dmso,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

nut1h <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin1h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut1h <- annotate_interactions(nut1h,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

intersect <- intersect_interactions(list(DMSO=dmso,Nutlin1h=nut1h))
gained_ints <- intersect$intersections$Nutlin1h
maintained_ints <- unique(intersect$intersections$`DMSO:Nutlin1h`)
lost_ints <- intersect$intersections$DMSO

source("../../scripts/integration_function.R")
link_gained <- p53_genes_interactions(ints=gained_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_gained <- link_gained[!grepl("non-annotated|\\.",link_gained$ensg),]
link_gained <- link_gained[link_gained$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_gained$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_gained <- left_join(link_gained,x)
link_gained <- link_gained[!is.na(link_gained$log2FC),]
link_gained <- link_gained %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_gained <- unique(link_gained)
link_gained$dynamism_interaction <- "gained"

link_maint <- p53_genes_interactions(ints=maintained_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_maint <- link_maint[!grepl("non-annotated|\\.",link_maint$ensg),]
link_maint <- link_maint[link_maint$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_maint$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_maint <- left_join(link_maint,x)
link_maint <- link_maint[!is.na(link_maint$log2FC),]
link_maint <- link_maint %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_maint <- unique(link_maint)
link_maint$dynamism_interaction <- "maintained"

link_lost <- p53_genes_interactions(ints=lost_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_lost <- link_lost[!grepl("non-annotated|\\.",link_lost$ensg),]
link_lost <- link_lost[link_lost$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_lost$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin1h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_lost <- left_join(link_lost,x)
link_lost <- link_lost[!is.na(link_lost$log2FC),]
link_lost <- link_lost %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_lost <- unique(link_lost)
link_lost$dynamism_interaction <- "lost"

link_res <- rbind(link_gained,link_maint,link_lost)
link_res <- link_res %>% group_by(p53) %>% filter(log2FC==max(log2FC))

link_res$log2FC <- NULL
write_rds(link_res,"../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1hvsDMSO.Rds")


## Distal integration DYNAMISM Nutlin10h vs. DMSO - Figure 6 ----
## functional p53 ----
functional_p53 <- read_rds("../../data/ComplementaryData/functional_p53.Rds")
functional_p53 <- functional_p53[grepl("Nutlin10h",functional_p53$condition),]

## Define promoters ----
promoters <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds")

# Integrate p53 with CREs ----
## Functional p53 at promoter ----
H3K27ac <- read_rds("../../data/ComplementaryData/H3K27ac/H3K27ac_Nutlin10h_merge_ChIP_input_peaks.narrowPeak")
active_promoters <- subsetByOverlaps(promoters,H3K27ac)
functional_p53_p <- subsetByOverlaps(functional_p53,active_promoters)

## Functional p53 at enhancer ----
functional_p53_e <- functional_p53[!functional_p53 %in% functional_p53_p,]
enhancers <- read_rds(file = "../../data/ComplementaryData/enhancers.Rds")
enhancers <- enhancers[grepl("Nutlin10h",enhancers$id),]
functional_p53_e <- subsetByOverlaps(functional_p53_e,enhancers)


DE_H3K27ac_Nutlin1h <-  read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")
DE_expr_Nutlin1h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin1h.Rds")
DE_expr_Nutlin10h <- read_rds("../../data/ComplementaryData/expression/DE_Expr_Nutlin10h.Rds")

## PCHi-C dynamism ----
dmso <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
dmso <- annotate_interactions(dmso,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

nut10h <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin10h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut10h <- annotate_interactions(nut10h,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

intersect <- intersect_interactions(list(DMSO=dmso,Nutlin10h=nut10h))
gained_ints <- intersect$intersections$Nutlin10h
maintained_ints <- unique(intersect$intersections$`DMSO:Nutlin10h`)
lost_ints <- intersect$intersections$DMSO

source("../../scripts/integration_function.R")
link_gained <- p53_genes_interactions(ints=gained_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_gained <- link_gained[!grepl("non-annotated|\\.",link_gained$ensg),]
link_gained <- link_gained[link_gained$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_gained$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_gained <- left_join(link_gained,x)
link_gained <- link_gained[!is.na(link_gained$log2FC),]
link_gained <- link_gained %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_gained <- unique(link_gained)
link_gained$dynamism_interaction <- "gained"


link_maint <- p53_genes_interactions(ints=maintained_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_maint <- link_maint[!grepl("non-annotated|\\.",link_maint$ensg),]
link_maint <- link_maint[link_maint$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_maint$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_maint <- left_join(link_maint,x)
link_maint <- link_maint[!is.na(link_maint$log2FC),]
link_maint <- link_maint %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_maint <- unique(link_maint)
link_maint$dynamism_interaction <- "maintained"

link_lost <- p53_genes_interactions(ints=lost_ints,p53.chip =functional_p53_e, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
link_lost <- link_lost[!grepl("non-annotated|\\.",link_lost$ensg),]
link_lost <- link_lost[link_lost$ensg %in% promoters$ensembl_gene_id,]
prom <- promoters[promoters$ensembl_gene_id %in% link_lost$ensg,]
x <- mergeByOverlaps(prom,DE_H3K27ac_Nutlin10h)
x <- data.frame(ensg=x$ensembl_gene_id,log2FC=x$log2FoldChange)
link_lost <- left_join(link_lost,x)
link_lost <- link_lost[!is.na(link_lost$log2FC),]
link_lost <- link_lost %>% group_by(p53) %>% filter(log2FC==max(log2FC))
link_lost <- unique(link_lost)
link_lost$dynamism_interaction <- "lost"

link_res <- rbind(link_gained,link_maint,link_lost)
link_res <- link_res %>% group_by(p53) %>% filter(log2FC==max(log2FC))

link_res$log2FC <- NULL
write_rds(link_res,"../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10hvsDMSO.Rds")
