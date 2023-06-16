# Figure 5 ----
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(HiCaptuRe)
library(GenomicRanges)
library(DESeq2)
library(VennDiagram)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

## Figure 5A ----
# Venn diagram distal target genes

nut1h <- read_rds("../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1h.Rds") # Tables resulting from distal_integration_PCHiC.R
nut10h <- read_rds("../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10h.Rds")

# visualize
grid.newpage()
draw.pairwise.venn(area1 = length(unique(nut1h$ensg)),    # Draw pairwise venn diagram
                   area2 = length(unique(nut10h$ensg)),
                   cross.area = sum(unique(nut1h$ensg) %in% unique(nut10h$ensg)),
                   fill = c("#37b578ff","#43377fff"),
                   category = c("Nut.1h","Nut.10h"),
                   rotation.degree = 180,cex = 2,
                   cat.pos=0)

## add how many genes are connected to 253 p53 from nut1h and nut10h interactomes in dmso 
gr1 <- makeGRangesFromDataFrame(nut1h %>% separate(p53,c("chr","start","end"),sep=":|-"))
gr2 <- makeGRangesFromDataFrame(nut10h %>% separate(p53,c("chr","start","end"),sep=":|-"))
p53_gr <- unique(c(gr1,gr2))

## interactions
dmso <- HiCaptuRe::load_interactions("../data/ComplementaryData/PCHiC/DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed")
dmso <- HiCaptuRe::annotate_interactions(dmso,annotation = "../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")

## promoters
promoters <- read_rds("../data/ComplementaryData/promoters_PCHiC.Rds")

## integrate
source("../scripts/integration_function.R")
link_dmso <- p53_genes_interactions(ints=dmso,p53.chip =p53_gr, promoter.annot = promoters,at.promoter=F,interaction.type = "P_OE")
upset(fromList(list(DMSO=unique(link_dmso$ensg),Nutlin1h=unique(nut1h$ensg),Nutlin10h=unique(nut10h$ensg))))



## Figure 5C ----
# Genomic distance distribution nearest genes vs. target genes of p53-bound enhancers

## nutlin1h 
tbl <- unique(nut1h[,colnames(nut1h) %in% c("ensg","p53")])

# promoters of genes of interest
promoters_genes <- promoters[promoters$ensembl_gene_id %in% tbl$ensg,]
tbl <- tbl[(tbl$ensg %in% promoters_genes$ensembl_gene_id),]
promoters_genes <- promoters_genes[match(tbl$ensg,promoters_genes$ensembl_gene_id),]

# p53 chip
p53_chip <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(tbl$p53,":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")

# compute distance distal target genes
gi <- GenomicInteractions(promoters_genes,p53_chip)
dist_p53_targetgene <- calculateDistances(gi)
target_genes <- gi@elementMetadata$anchor1.ensembl_gene_id
dist_p53_targetgene <- dist_p53_targetgene[!is.na(dist_p53_targetgene)]

# compute distal nearest target gene
promoters_all <- read_rds("../data/ComplementaryData/promoters_GRCh37.87.Rds")
nearest_p53chip <- nearest(x = unique(p53_chip), subject = promoters_all, ignore.strand=T)
gi2 <- GenomicInteractions(promoters_all[nearest_p53chip],unique(p53_chip))
nearest_genes <- gi2@elementMetadata$anchor1.ensembl_gene_id
dist_p53_nearestgene <- calculateDistances(gi2)

to.plot <- data.frame(variable=c(rep("target_gene",length(dist_p53_targetgene)),rep("nearest_gene",length(dist_p53_nearestgene))),
                      value=c(dist_p53_targetgene,dist_p53_nearestgene))
to.plot$variable<-factor(to.plot$variable,levels = c("nearest_gene","target_gene"))
median.to.plot <- to.plot %>% group_by(variable) %>% summarize(med_value=median(value,na.rm=T))

p1 <- ggplot(to.plot,aes(log(value),col=variable)) + 
  geom_density(aes(fill=variable),alpha=0.3) +
  scale_color_manual(values=c("darkgrey","#37b578ff"),labels=c(paste0("nearest gene (",length(unique(nearest_genes)),")"),paste0("target gene (",length(unique(target_genes)),")"))) +
  scale_fill_manual(values=c("darkgrey","#37b578ff"),labels=c(paste0("nearest gene (",length(unique(nearest_genes)),")"),paste0("target gene (",length(unique(target_genes)),")"))) +
  labs(x="Distance log(bp)",y="Density",fill="",color="",title = "Nutlin1h interactions") +
  theme_classic2() + 
  geom_vline(xintercept = median(log(to.plot[to.plot$variable=="target_gene",]$value)),color="#37b578ff") +
  geom_vline(xintercept = median(log(to.plot[to.plot$variable=="nearest_gene",]$value)),color="darkgrey") +
  geom_text(data=median.to.plot,mapping = aes(x=log(med_value),y=0.2,label=paste0("median=",round(med_value,3)),col=variable),color="black",angle=270) +
  annotate(geom="text",x=18,y=0.4,label=c(paste0("mean target = ",round(mean(to.plot[to.plot$variable=="target_gene",]$value))))) +
  annotate(geom="text",x=18,y=0.35,label=c(paste0("mean nearest = ",round(mean(to.plot[to.plot$variable=="nearest_gene",]$value))))) + 
  xlim(0,20)

## nutlin10h 
tbl <- unique(nut10h[,colnames(nut10h) %in% c("ensg","p53")])

# promoters of genes of interest
promoters_genes <- promoters[promoters$ensembl_gene_id %in% tbl$ensg,]
tbl <- tbl[(tbl$ensg %in% promoters_genes$ensembl_gene_id),]
promoters_genes <- promoters_genes[match(tbl$ensg,promoters_genes$ensembl_gene_id),]

# p53 chip
p53_chip <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(tbl$p53,":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")

# compute distance distal target genes
gi <- GenomicInteractions(promoters_genes,p53_chip)
dist_p53_targetgene <- calculateDistances(gi)
target_genes <- gi@elementMetadata$anchor1.ensembl_gene_id
dist_p53_targetgene <- dist_p53_targetgene[!is.na(dist_p53_targetgene)]

# compute distance nearest target genes
promoters_all <- read_rds("../data/ComplementaryData/promoters_GRCh37.87.Rds")
nearest_p53chip <- nearest(x = unique(p53_chip), subject = promoters_all, ignore.strand=T)
gi2 <- GenomicInteractions(promoters_all[nearest_p53chip],unique(p53_chip))
nearest_genes <- gi2@elementMetadata$anchor1.ensembl_gene_id
dist_p53_nearestgene <- calculateDistances(gi2)

to.plot <- data.frame(variable=c(rep("target_gene",length(dist_p53_targetgene)),rep("nearest_gene",length(dist_p53_nearestgene))),
                      value=c(dist_p53_targetgene,dist_p53_nearestgene))
to.plot$variable<-factor(to.plot$variable,levels = c("nearest_gene","target_gene"))
median.to.plot <- to.plot %>% group_by(variable) %>% summarize(med_value=median(value,na.rm=T))

p2 <- ggplot(to.plot,aes(log(value),col=variable)) + 
  geom_density(aes(fill=variable),alpha=0.3) +
  scale_color_manual(values=c("darkgrey","#43377fff"),labels=c(paste0("nearest gene (",length(unique(nearest_genes)),")"),paste0("target gene (",length(unique(target_genes)),")"))) +
  scale_fill_manual(values=c("darkgrey","#43377fff"),labels=c(paste0("nearest gene (",length(unique(nearest_genes)),")"),paste0("target gene (",length(unique(target_genes)),")"))) +
  labs(x="Distance log(bp)",y="Density",fill="",color="",title = "Nutlin10h interactions") +
  theme_classic2() + 
  geom_vline(xintercept = median(log(to.plot[to.plot$variable=="target_gene",]$value)),color="#43377fff") +
  geom_vline(xintercept = median(log(to.plot[to.plot$variable=="nearest_gene",]$value)),color="darkgrey") +
  geom_text(data=median.to.plot,mapping = aes(x=log(med_value),y=0.2,label=paste0("median=",round(med_value,3)),col=variable),color="black",angle=270) +
  annotate(geom="text",x=18,y=0.4,label=c(paste0("mean target = ",round(mean(to.plot[to.plot$variable=="target_gene",]$value))))) +
  annotate(geom="text",x=18,y=0.35,label=c(paste0("mean nearest = ",round(mean(to.plot[to.plot$variable=="nearest_gene",]$value))))) + 
  xlim(0,20)

#visualize
p1
p2



## Figure 5D ----
# GSEA bubble plot p53 target genes


nut1h <- read_rds("../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1h.Rds") # Tables resulting from distal_integration_PCHiC.R
nut10h <- read_rds("../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin10h.Rds")
nut1h_p <- read_rds("../data/ComplementaryData/PCHiC_integration/proximal_integration_Nutlin1h.Rds") # Tables from proximal_integration_PCHiC.R
nut10h_p <- read_rds("../data/ComplementaryData/PCHiC_integration/proximal_integration_Nutlin10h.Rds")

## GFF promoters
promoters_all <- read_rds("../data/ComplementaryData/promoters_GRCh37.87.Rds")

## PCHiC promoters
promoters <- read_rds("../data/ComplementaryData/promoters_PCHiC.Rds")

## return target genes and nearest genes for 2 conditions
list_link <- list(nutlin1h=nut1h,nutlin10h=nut10h)
list_genes <- lapply(list_link,function(tbl){
  
  tbl <- unique(tbl[,colnames(tbl) %in% c("ensg","p53")])
  
  # target gene
  target_genes <- tbl$ensg
  
  # nearest gene
  p53_chip <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(tbl$p53,":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
  nearest_p53chip <- nearest(x = unique(p53_chip), subject = promoters_all, ignore.strand=T)
  gi <- GenomicInteractions(promoters_all[nearest_p53chip],unique(p53_chip))
  nearest_genes <- gi@elementMetadata$anchor1.ensembl_gene_id
  
  return(list(target=target_genes,nearest=nearest_genes))
})

# proximal gene
list_genes$nutlin1h$promoter <- unique(nut1h_p$target_gene)
names(list_genes$nutlin1h) <- paste0("nutlin1h_",names(list_genes$nutlin1h))
list_genes$nutlin10h$promoter <- unique(nut10h_p$target_gene)
names(list_genes$nutlin10h) <- paste0("nutlin10h_",names(list_genes$nutlin10h))

to.plot <- c(list_genes$nutlin1h,list_genes$nutlin10h)
to.plot <- lapply(to.plot,unique)

# ensembl to gene id
genes_annotation <- as.data.frame(promoters@elementMetadata)

df <- data.frame(matrix(nrow = 0, ncol = 0)) 
categories <- list(H="H",C1="C1",C2="C2",C3="C3",C4="C4",C5="C5",C6="C6",C7="C7",C8="C8")
res <- lapply(categories,function(cat){
  print(cat)
  for(i in 1:length(to.plot)){
    gene.set <- msigdbr(species = "Homo sapiens",category = cat)
    msigdbr_t2g <- gene.set %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
    
    res.gsea <- enricher(gene =unlist(to.plot[[i]]), TERM2GENE = msigdbr_t2g,universe =promoters_all$ensembl_gene_id )
    res.gsea <- as.data.frame(res.gsea)
    
    aux <- data.frame(matrix(nrow = 0, ncol = 0)) 
    if(nrow(res.gsea)!=0){
      res.gsea$type <- names(to.plot)[i]
      df <- rbind(df,res.gsea)
    }
    print(names(to.plot)[i])
    print(df)
  }
  return(df)
})
for(i in 1:length(res)){
  if(nrow(res[[i]]!=0)){
    res[[i]]$category <- names(res)[i]
  }
}
df <- do.call("rbind",res)

df <- df[df$category %in% c("C2","C5","C6"),]
x <- df %>% group_by(type) %>% arrange(p.adjust) %>% slice(1:5)
plot_gsea <- df[df$ID %in% unique(x$ID),]
plot_gsea <-plot_gsea %>% arrange(type)
plot_gsea <- plot_gsea %>% separate(type,c("condition","group"),sep="_",remove = F)
plot_gsea$condition <- factor(plot_gsea$condition,levels=c("nutlin1h","nutlin10h"))
plot_gsea$group <- factor(plot_gsea$group,levels=c("nearest","promoter","target"))
plot_gsea <- plot_gsea %>% arrange(desc(condition),group)
plot_gsea$ID <- factor(plot_gsea$ID,levels=unique(plot_gsea$ID))

# visualize
ggplot(plot_gsea,aes(type,ID,size=Count,color=p.adjust)) + 
  geom_point() + 
  theme_bw() +
  facet_grid(cols=vars(condition),scales="free",space = "free_y") +
  labs(x = "p53 target gene groups",y="MSigDB gene set",size="#genes",color="p-value") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_colour_viridis(discrete = F,begin = 0.8)


## Figure 5E/H ----
# H3K27ac of CREs: enhancer-bound p53, promoter-bound p53 and control

DE_H3K27ac_Nutlin1h <- read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin1h.Rds")
DE_H3K27ac_Nutlin10h <- read_rds("../data/ComplementaryData/H3K27ac/DE_H3K27ac_Nutlin10h.Rds")

## H3K27ac log2FC at promoters nutlin1h
gr_p <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut1h_p$promoter_target_gene),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_p,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_p))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df1 <- unique(df)[,2:3]

## H3K27ac log2FC at enhancers nutlin1h
gr_e <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut1h$enhancer),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_e,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_e))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df2 <- unique(df)[,2:3]

## H3K27ac log2FC background 
gr_bg <- subsetByOverlaps(DE_H3K27ac_Nutlin1h,c(gr_p,gr_e),invert=T)
df <- data.frame(functional_element=rep("background",length(gr_bg)),log2FC=gr_bg$log2FoldChange,element=as.character(gr_bg))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df3 <- unique(df)[,2:3]

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC))

p1 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0, linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape=NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("enhancer bound by p53 \n (n=",length(gr_e),")"),
                            paste0("promoter bound by p53 \n (n=",length(gr_p),")"),
                            paste0("background acetylation"))
  ) +
  labs(x="p53-bound regulatory element",y="mean H3K27ac log2FC",col="") +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,2)),nudge_y = 0.15) +
  ylim(-2,3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

## H3K27ac log2FC at promoters nutlin10h
gr_p <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut10h_p$promoter_target_gene),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_p,DE_H3K27ac_Nutlin10h)
df <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_p))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df1 <- unique(df)[,2:3]

## H3K27ac log2FC at enhancers nutlin10h
gr_e <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut10h$enhancer),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_e,DE_H3K27ac_Nutlin10h)
df <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_e))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df2 <- unique(df)[,2:3]

## H3K27ac log2FC background 
gr_bg <- subsetByOverlaps(DE_H3K27ac_Nutlin10h,c(gr_p,gr_e),invert=T)
df <- data.frame(functional_element=rep("background",length(gr_bg)),log2FC=gr_bg$log2FoldChange,element=as.character(gr_bg))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df3 <- unique(df)[,2:3]

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC))

p2 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0, linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape=NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("enhancer bound by p53 \n (n=",length(gr_e),")"),
                            paste0("promoter bound by p53 \n (n=",length(gr_p),")"),
                            paste0("background acetylation")
  )
  ) +
  labs(x="p53-bound regulatory element",y="mean H3K27ac log2FC",col="") +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,2)),nudge_y = 0.15) +
  # ylim(-2,3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

#visualize
p1
p2


## Figure 5F/I ----
# H3K27ac of target promoter to enhancer-bound p53, promoter-bound p53 and control

## H3K27ac log2FC at promoter bound p53 target genes 
gr_p <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut1h_p$promoter_target_gene),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_p,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_p))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df1 <- unique(df)

## H3K27ac log2FC at enhancer bound p53 target genes 
gr_e <- promoters[promoters$ensembl_gene_id %in% nut1h$ensg,]
merge <- mergeByOverlaps(gr_e,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_e))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df2 <- unique(df)

## H3K27ac log2FC at background promoters 
gr_bg <- promoters[!promoters$ensembl_gene_id %in% c(nut1h$ensg,nut1h_p$target_gene),]
merge <- mergeByOverlaps(gr_bg,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("background",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_bg))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df3 <- unique(df)

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC))

p1 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape = NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("target to enhancer-p53 \n (n=",length(unique(df2$element)),")"),
                            paste0("target to promoter-p53 \n (n=",length(unique(df1$element)),")"),
                            paste0("all non-target \n (n=",length(unique(df3$element)),")")
  )
  ) +
  labs(x="target genes to p53",y="mean H3K27ac log2FC at target gene promoters",col="")  + 
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,2)),nudge_y = 0.15) +
  # ylim(-2,4) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 


## H3K27ac log2FC at promoter bound p53 target genes nut1h
gr_p <- makeGRangesFromDataFrame(as.data.frame(do.call("rbind",str_split(unique(nut10h_p$promoter_target_gene),":|-"))),seqnames.field = "V1",start.field = "V2",end.field = "V3")
merge <- mergeByOverlaps(gr_p,DE_H3K27ac_Nutlin10h)
df <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_p))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df1 <- unique(df)

## H3K27ac log2FC at enhancer bound p53 target genes
gr_e <- promoters[promoters$ensembl_gene_id %in% nut10h$ensg,]
merge <- mergeByOverlaps(gr_e,DE_H3K27ac_Nutlin1h)
df <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_e))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df2 <- unique(df)

## H3K27ac log2FC at background promoters
gr_bg <- promoters[!promoters$ensembl_gene_id %in% c(nut10h$ensg,nut10h_p$target_gene),]
merge <- mergeByOverlaps(gr_bg,DE_H3K27ac_Nutlin10h)
df <- data.frame(functional_element=rep("background",nrow(merge)),log2FC=merge$log2FoldChange,element=as.character(merge$gr_bg))
df <- df %>% group_by(element) %>% summarize(log2FC=mean(log2FC,rm.na=T),functional_element=functional_element)
df3 <- unique(df)

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC))

p2 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape = NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("target to enhancer-p53 \n (n=",length(unique(df2$element)),")"),
                            paste0("target to promoter-p53 \n (n=",length(unique(df1$element)),")"),
                            paste0("all non-target \n (n=",length(unique(df3$element)),")")
  )
  ) +
  labs(x="target genes to p53",y="mean H3K27ac log2FC at target gene promoters",col="")  + 
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,2)),nudge_y = 0.15) +
  # ylim(-2,4) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

#visualize
p1
p2


## Figure 5G/J ----
# Expression of target genes to enhancer-bound p53, promoter-bound p53 and control

expr_1h <- read_rds("../data/DE_Expr_Nutlin1h.Rds")
expr_10h <- read_rds("../data/DE_Expr_Nutlin10h.Rds")

# p53 at promoter
merge <- expr_1h[expr_1h$ensg %in% nut1h_p$target_gene,]
df1 <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange)

# p53 at enhancer
merge <- expr_1h[expr_1h$ensg %in% nut1h$ensg,]
df2 <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange)

merge <- expr_1h[expr_1h$ensg %in% promoters$ensembl_gene_id[!promoters$ensembl_gene_id %in% c(nut1h$ensg,nut1h_p$target_gene)],]
df3 <- data.frame(functional_element=rep("background",nrow(merge)),log2FC=merge$log2FoldChange)

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC,na.rm=T))

p1 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape = NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#37b578ff","#37b578ff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("target to enhancer-p53 \n (n=",nrow(df2),")"),
                            paste0("target to promoter-p53 \n (n=",nrow(df1),")"),
                            paste0("all non-target \n (n=",nrow(df3),")")
  )
  ) +
  labs(x="target genes to p53",y="Expression 1h vs. 0h log2FC",col="") +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,3)),nudge_y = 0.3) +
  # ylim(-3,3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

# p53 at promoter
merge <- expr_10h[expr_10h$ensg %in% nut10h_p$target_gene,]
df1 <- data.frame(functional_element=rep("promoter",nrow(merge)),log2FC=merge$log2FoldChange)

# p53 at enhancer
merge <- expr_10h[expr_10h$ensg %in% nut10h$ensg,]
df2 <- data.frame(functional_element=rep("enhancer",nrow(merge)),log2FC=merge$log2FoldChange)

merge <- expr_10h[expr_10h$ensg %in% promoters$ensembl_gene_id[!promoters$ensembl_gene_id %in% c(nut10h$ensg,nut10h_p$target_gene)],]
df3 <- data.frame(functional_element=rep("background",nrow(merge)),log2FC=merge$log2FoldChange)

to.plot <- rbind(df1,df2,df3)
to.plot$functional_element <- factor(to.plot$functional_element,levels=c("enhancer","promoter","background"))

my_comparisons <- list(c("enhancer","background"),c("promoter","background"))
to.plot.annot <- to.plot %>% group_by(functional_element) %>% summarise(med=median(log2FC,na.rm=T))

p2 <- ggplot(to.plot,aes(functional_element,log2FC)) + 
  geom_violin(aes(col=functional_element),width=0.6) +
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey") +
  geom_boxplot(aes(fill=functional_element),width=0.2,outlier.shape = NA) + 
  theme_classic2() + 
  scale_colour_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("#43377fff","#43377fff","darkgrey"),guide=NULL) +
  scale_x_discrete(labels=c(paste0("target to enhancer-p53 \n (n=",nrow(df2),")"),
                            paste0("target to promoter-p53 \n (n=",nrow(df1),")"),
                            paste0("all non-target \n (n=",nrow(df3),")")
  )
  ) +
  labs(x="target genes to p53",y="Expression 10h vs. 0h log2FC",col="") +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test") +
  geom_text(data = to.plot.annot,aes(functional_element,med,label=round(med,3)),nudge_y = 0.3) +
  ylim(-3,3) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

#visualise
p1
p2


## Figure 5K ----
# qRT-PCR relative expression of distal target genes

x <- read_rds("../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression.Rds")


x <- x %>% group_by(Cond,variable) %>% summarize(mean=mean(value),sd=sd(value))
x <- x[grepl("WT",x$Cond),]
x$Cond <- factor(x$Cond,levels = c("WTNut0h","WTNut1h","WTNut10h"))
x <- x[!x$variable %in% c("PPP1CB","CHD4","ACTR8","TP53","TGFA","BAX","CDKN1A"),]
p<- ggplot(x, aes(x=variable, y=mean, fill=Cond)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  ggpubr::theme_classic2() + 
  scale_fill_manual(values = c("yellow","green","purple"))

# visualise
p

### T-test between timepoints within genes
## WTNutlin10h vs. WTDMSO
df <- as.data.frame(matrix(nrow=4))[,-1]
for(gene in unique(x$variable)){
  print(gene)
  res <- t.test(x[x$Cond=="WTNut0h"&x$variable==gene,]$value,x[x$Cond=="WTNut10h"&x$variable==gene,]$value,alternative = "two.sided")
  testres <- res[[1]]
  dfres <- res[[2]]
  pvalue <- res[[3]]
  vec <- unlist(c(gene,testres,dfres,pvalue))
  df <- rbind(df,as.data.frame(t(as.data.frame(vec))))
}
rownames(df) <- 1:nrow(df)
colnames(df) <- c("gene","t.test.stat","t.test.df","t.test.pval")
df$comparison <- "WTNut10h.WTNut0h"
df1 <- df

## WTNutlin1h vs. WTDMSO
df <- as.data.frame(matrix(nrow=4))[,-1]
for(gene in unique(x$variable)){
  print(gene)
  res <- t.test(x[x$Cond=="WTNut0h"&x$variable==gene,]$value,x[x$Cond=="WTNut1h"&x$variable==gene,]$value,alternative = "two.sided")
  testres <- res[[1]]
  dfres <- res[[2]]
  pvalue <- res[[3]]
  vec <- unlist(c(gene,testres,dfres,pvalue))
  df <- rbind(df,as.data.frame(t(as.data.frame(vec))))
}
rownames(df) <- 1:nrow(df)
colnames(df) <- c("gene","t.test.stat","t.test.df","t.test.pval")
df$comparison <- "WTNut1h.WTNut0h"
df2 <- df
