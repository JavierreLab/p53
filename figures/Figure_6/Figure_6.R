# Figure 6 ----
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(HiCaptuRe)
library(GenomicRanges)
library(DESeq2)
library(GenomicInteractions)


## Figure 6A ----

# promoters
promoters_all <- read_rds("../../data/ComplementaryData/promoters_GRCh37.87.Rds") # promoters -1000/+200 all genomic GTF GRCh37v87
promoters <- read_rds("../../data/ComplementaryData/promoters_PCHiC.Rds") # promoters computes as -1000/+200 from PCHiC annotation

# distal integration data
nut1h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1hvsDMSO.Rds") # Tables resulting from distal_integration_PCHiC.R
nut10h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10hvsDMSO.Rds")

# nutlin1h
df <- as.data.frame(table(unique(nut1h[,colnames(nut1h) %in% c("dynamism_interaction","bait","oe")])$dynamism_interaction))
df$type <- "p53"

dmso_ints <- load_interactions(file = "../../data/ComplementaryData/PCHiC/DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed")
dmso_ints <- annotate_interactions(dmso_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")
dmso_ints_vec <- paste0(as.character(anchorOne(dmso_ints)),as.character(anchorTwo(dmso_ints)))

nut1h_ints <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin1h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut1h_ints <- annotate_interactions(nut1h_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")
nut1h_ints_vec <- paste0(as.character(anchorOne(nut1h_ints)),as.character(anchorTwo(nut1h_ints)))

as.vector(table(nut1h_ints_vec %in% dmso_ints_vec))
df2 <- data.frame(Var1=c("gained","maintained","lost"),Freq=c(as.vector(table(nut1h_ints_vec %in% dmso_ints_vec)),as.vector(table(dmso_ints_vec %in% nut1h_ints_vec))[1]),type="total")

dmso_poe <- dmso_ints[dmso_ints$int=="P_OE"]
nut1h_poe <- nut1h_ints[nut1h_ints$int=="P_OE"]
dmso_ints_poe <- paste0(as.character(anchorOne(dmso_poe)),as.character(anchorTwo(dmso_poe)))
nut1h_ints_poe <- paste0(as.character(anchorOne(nut1h_poe)),as.character(anchorTwo(nut1h_poe)))

df3 <- data.frame(Var1=c("gained","maintained","lost"),Freq=c(as.vector(table(nut1h_ints_poe %in% dmso_ints_poe)),as.vector(table(dmso_ints_poe %in% nut1h_ints_poe))[1]),type="P_OE")

to.plot <- rbind(df,df2,df3)
to.plot$type <- factor(to.plot$type,levels = c("total","P_OE","p53"))
to.plot$Var1 <- factor(as.character(to.plot$Var1),levels = c("maintained","gained","lost"))

p1 <- ggplot(to.plot,aes(type,Freq,fill=Var1)) +
  geom_col(position= "fill") +
  labs(x = "Interaction set", y = "Interactions (%)", fill = NULL) +
  theme_classic2() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(~type,scales = "free_x") +
  scale_fill_manual(values=c("darkgrey","grey","lightgrey")) +
  geom_text(aes(label=Freq),position=position_fill(vjust = 0.5))

# nutlin10h
df <- as.data.frame(table(unique(nut10h[,colnames(nut10h) %in% c("dynamism_interaction","bait","oe")])$dynamism_interaction))
df$type <- "p53"

dmso_ints <- load_interactions(file = "../../data/ComplementaryData/PCHiC/DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed")
dmso_ints <- annotate_interactions(dmso_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")
dmso_ints_vec <- paste0(as.character(anchorOne(dmso_ints)),as.character(anchorTwo(dmso_ints)))

nut10h_ints <- load_interactions(file = paste0("../../data/ComplementaryData/PCHiC/Nutlin10h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed"))
nut10h_ints <- annotate_interactions(nut10h_ints,annotation = "../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")
nut10h_ints_vec <- paste0(as.character(anchorOne(nut10h_ints)),as.character(anchorTwo(nut10h_ints)))

as.vector(table(nut10h_ints_vec %in% dmso_ints_vec))
df2 <- data.frame(Var1=c("gained","maintained","lost"),Freq=c(as.vector(table(nut10h_ints_vec %in% dmso_ints_vec)),as.vector(table(dmso_ints_vec %in% nut10h_ints_vec))[1]),type="total")

dmso_poe <- dmso_ints[dmso_ints$int=="P_OE"]
nut10h_poe <- nut10h_ints[nut10h_ints$int=="P_OE"]
dmso_ints_poe <- paste0(as.character(anchorOne(dmso_poe)),as.character(anchorTwo(dmso_poe)))
nut10h_ints_poe <- paste0(as.character(anchorOne(nut10h_poe)),as.character(anchorTwo(nut10h_poe)))

df3 <- data.frame(Var1=c("gained","maintained","lost"),Freq=c(as.vector(table(nut10h_ints_poe %in% dmso_ints_poe)),as.vector(table(dmso_ints_poe %in% nut10h_ints_poe))[1]),type="P_OE")

to.plot <- rbind(df,df2,df3)
to.plot$type <- factor(to.plot$type,levels = c("total","P_OE","p53"))
to.plot$Var1 <- factor(as.character(to.plot$Var1),levels = c("maintained","gained","lost"))

p2 <- ggplot(to.plot,aes(type,Freq,fill=Var1)) +
  geom_col(position= "fill") +
  labs(x = "Interaction set", y = "Interactions (%)", fill = NULL) +
  theme_classic2() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(~type,scales = "free_x") +
  scale_fill_manual(values=c("darkgrey","grey","lightgrey")) +
  geom_text(aes(label=Freq),position=position_fill(vjust = 0.5))

# visualise
p1
p2


## Figure 6B ----

# distal integration data
nut1h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1hvsDMSO.Rds") # Tables resulting from distal_integration_PCHiC.R
nut10h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10hvsDMSO.Rds")

list_genes <- list(gained_interactions=nut1h[nut1h$dynamism_interaction=="gained",]$ensg,
                   maintained_interactions=nut1h[nut1h$dynamism_interaction=="maintained",]$ensg,
                   lost_interactions=nut1h[nut1h$dynamism_interaction=="lost",]$ensg)
p1 <- UpSetR::upset(UpSetR::fromList(list_genes),mainbar.y.label = "distal target genes to p53")

list_genes <- list(gained_interactions=nut10h[nut10h$dynamism_interaction=="gained",]$ensg,
                   maintained_interactions=nut10h[nut10h$dynamism_interaction=="maintained",]$ensg,
                   lost_interactions=nut10h[nut10h$dynamism_interaction=="lost",]$ensg)
p2 <- UpSetR::upset(UpSetR::fromList(list_genes),mainbar.y.label = "distal target genes to p53")

p1
p2


list_p53 <- list(gained_interactions=nut1h[nut1h$dynamism_interaction=="gained",]$p53,
                 maintained_interactions=nut1h[nut1h$dynamism_interaction=="maintained",]$p53,
                 lost_interactions=nut1h[nut1h$dynamism_interaction=="lost",]$p53)
p3 <- UpSetR::upset(UpSetR::fromList(list_p53),mainbar.y.label = "p53")

list_p53 <- list(gained_interactions=nut10h[nut10h$dynamism_interaction=="gained",]$p53,
                 maintained_interactions=nut10h[nut10h$dynamism_interaction=="maintained",]$p53,
                 lost_interactions=nut10h[nut10h$dynamism_interaction=="lost",]$p53)
p4 <- UpSetR::upset(UpSetR::fromList(list_p53),mainbar.y.label = "p53")

p3
p4

## Figure 6C ----

# distal integration data
nut1h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1hvsDMSO.Rds") # Tables resulting from distal_integration_PCHiC.R
nut10h <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10hvsDMSO.Rds")

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

##nutlin1h
gained.genes <- unique(nut1h[nut1h$dynamism_interaction=="gained",]$ensg)
maintained.genes <- unique(nut1h[nut1h$dynamism_interaction=="maintained",]$ensg)
lost.genes <- unique(nut1h[nut1h$dynamism_interaction=="lost",]$ensg)

# universe
universe <- promoters$ensembl_gene_id

list_genes <- list(gained=gained.genes, maintained=maintained.genes,
                   lost=lost.genes)
df <- data.frame(matrix(nrow = 0, ncol = 0)) 
for(cat in c("H","C1","C2","C3","C4","C5","C6","C7","C8")){
  print(cat)
  gene.set <- msigdbr(species = "Homo sapiens",category = cat)
  msigdbr_t2g <- gene.set %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
  aux <- data.frame(matrix(nrow = 0, ncol = 0)) 
  for(i in 1:length(list_genes)){
    print(names(list_genes)[i])
    res <- as.data.frame(enricher(gene =unlist(list_genes[i]), TERM2GENE = msigdbr_t2g,universe = universe))
    if(nrow(res)!=0){
      res$group <- names(list_genes)[i]
      res$category <- cat
      aux <- rbind(aux,res)
    }
  }
  df <- rbind(df,aux)
}
df <- df[!df$category %in% c("C1"),]
x <- df %>% group_by(group) %>% arrange(p.adjust) %>% slice(1:5)
df <- df[df$ID %in% unique(x$ID),]
df$group <- factor(df$group,levels=c("maintained","gained","lost"))
df <- df %>% arrange(group,ID)
df$ID <- factor(df$ID,levels=unique(df$ID))

# visualise
ggplot(df,aes(group,ID,size=Count,color=p.adjust)) + 
  geom_point() + 
  theme_classic2() +
  viridis::scale_colour_viridis(discrete = F,begin = 0.8) +
  facet_grid(cols=vars(group),scales="free",space = "free_y") +
  labs(x = "",y="Gene Set",size="# genes",color="p-value")

## nutlin10h
gained.genes <- unique(nut10h[nut10h$dynamism_interaction=="gained",]$ensg)
maintained.genes <- unique(nut10h[nut10h$dynamism_interaction=="maintained",]$ensg)
lost.genes <- unique(nut10h[nut10h$dynamism_interaction=="lost",]$ensg)

# universe
universe <- promoters$ensembl_gene_id

list_genes <- list(gained=gained.genes, maintained=maintained.genes,
                   lost=lost.genes)
df <- data.frame(matrix(nrow = 0, ncol = 0)) 
for(cat in c("H","C1","C2","C3","C4","C5","C6","C7","C8")){
  print(cat)
  gene.set <- msigdbr(species = "Homo sapiens",category = cat)
  msigdbr_t2g <- gene.set %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
  aux <- data.frame(matrix(nrow = 0, ncol = 0)) 
  for(i in 1:length(list_genes)){
    print(names(list_genes)[i])
    res <- as.data.frame(enricher(gene =unlist(list_genes[i]), TERM2GENE = msigdbr_t2g,universe = universe))
    if(nrow(res)!=0){
      res$group <- names(list_genes)[i]
      res$category <- cat
      aux <- rbind(aux,res)
    }
  }
  df <- rbind(df,aux)
}
df <- df[!df$category %in% c("C1"),]
x <- df %>% group_by(group) %>% arrange(p.adjust) %>% slice(1:5)
df <- df[df$ID %in% unique(x$ID),]
df$group <- factor(df$group,levels=c("maintained","gained","lost"))
df <- df %>% arrange(group,ID)
df$ID <- factor(df$ID,levels=unique(df$ID))

# visualise
ggplot(df,aes(group,ID,size=Count,color=p.adjust)) + 
  geom_point() + 
  theme_classic2() +
  viridis::scale_colour_viridis(discrete = F,begin = 0.8) +
  facet_grid(cols=vars(group),scales="free",space = "free_y") +
  labs(x = "",y="Gene Set",size="# genes",color="p-value")



