# Figure 4 ----
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(DESeq2)
library(HiCaptuRe)
library(VennDiagram) 
library(corrplot)
library(data.table)


## Figure 4A ----
pkm <- load_interactions(file = "../../data/ComplementaryData/PCHiC/peakmatrix_recalibrated2WTDMSO.txt") # chicago output
pkm2 <- as.data.frame(pkm@elementMetadata[,6:c(ncol(pkm@elementMetadata)-1)])

cor_matrix <- cor(pkm2)
colnames(cor_matrix)[grepl("Nutlin3a_1h",colnames(cor_matrix))] <- gsub("CS_Nutlin3a_1h_BR","Nut.1h",colnames(cor_matrix)[grepl("Nutlin3a_1h",colnames(cor_matrix))])
colnames(cor_matrix)[grepl("Nutlin3a_10h",colnames(cor_matrix))] <- gsub("CS_Nutlin3a_10h_BR","Nut.10h",colnames(cor_matrix)[grepl("Nutlin3a_10h",colnames(cor_matrix))])
colnames(cor_matrix)[grepl("DMSO",colnames(cor_matrix))] <- gsub("CS_DMSO_BR","Nut.0h",colnames(cor_matrix)[grepl("DMSO",colnames(cor_matrix))])
rownames(cor_matrix) <- colnames(cor_matrix)

col <- circlize::colorRamp2(breaks = c(0.2,0.5,1),colors = c("#ffffff","#dee0eb","#23447f"))

# visualize
ComplexHeatmap::Heatmap(cor_matrix,col=col,heatmap_legend_param = list(title = "Pearson Cor."))


## Figure 4B ----
pkm <- load_interactions(file = "../../data/ComplementaryData/PCHiC/peakmatrix_recalibrated2WTDMSO.txt") # chicago output
pkm2 <- as.data.frame(pkm@elementMetadata[,6:c(ncol(pkm@elementMetadata)-1)])

pca_res <- prcomp(t(pkm2))
res <- as.data.frame(pca_res$x)
res$sample <- rownames(res)
res$condition <- gsub("_BR[1-2]","",res$sample)
# visualize
ggplot(res,aes(PC1,PC2,col=condition))+ 
  geom_point(size=5)+
  # ggforce::geom_mark_ellipse(aes(fill=condition)) +
  theme_classic2() +
  scale_color_manual(values=c("#fde725ff","#43377fff","#37b578ff")) +
  scale_fill_manual(values=c("#fde725ff","#43377fff","#37b578ff")) +
  xlab(paste0("PC1 (",summary(pca_res)$importance[2,1]*100,"%)")) +
  ylab(paste0("PC2 (",summary(pca_res)$importance[2,2]*100,"%)")) + 
  labs(fill="",color="")


## Figure 4C ----

pkm <- load_interactions(file = "../../data/ComplementaryData/PCHiC/peakmatrix_merged_recalibrated2WTDMSO.txt") # chicago output
pkm <- as.data.frame(pkm@elementMetadata[,grepl("CS_",colnames(pkm@elementMetadata))])
pkm$id <- 1:nrow(pkm)

dmso <- pkm[pkm$CS_DMSO>5,]$id
nut1h <- pkm[pkm$CS_Nutlin3a_1h>5,]$id
nut10h <- pkm[pkm$CS_Nutlin3a_10h>5,]$id


# visualize
grid.newpage()
draw.pairwise.venn(area1 = length(unique(dmso)),
                   area2 = length(unique(nut1h)),
                   cross.area = sum(dmso %in% nut1h),
                   fill = c("#fde725ff","#37b578ff"),
                   category = c("Nut.0h","Nut.1h"),
                   cat.pos=0)
grid.newpage()
draw.pairwise.venn(area1 = length(unique(dmso)),
                   area2 = length(unique(nut10h)),
                   cross.area = sum(dmso %in% nut10h),
                   fill = c("#fde725ff","#43377fff"),
                   rotation.degree = 180,
                   category = c("Nut.0h","Nut.10h"),
                   cat.pos = 0)

## Figure 4E ----

wd <- "../../data/ComplementaryData/PCHiC/"

clean_ibeds <- list.files(path = paste0(wd), pattern = "recalibrated2WTDMSO_cutoff_5.ibed",full.names = T) # chicago output

interactions_list <- list()
PP=NULL
POE=NULL
short_distances= NULL

for (i in clean_ibeds)
{
  name <- basename(i)
  
  a <- HiCaptuRe::load_interactions(i)
  a <- HiCaptuRe::annotate_interactions(a, annotation ="../../data/ComplementaryData/PCHiC/baits_coordinates_annotated_ensembl_gene_id_GRCh37_87.bed")
  a_cis <- a[is.cis(a)]
  
  a_distances <- HiCaptuRe::distance_summary(interactions = a_cis, sample = name)
  
  a_short_distances <- a_distances
  short_distances=rbind(short_distances,a_short_distances)
  
  d_PP <- data.frame(a_short_distances[grep("P_P", a_short_distances$int), ])
  d_POE <- data.frame(a_short_distances[grep("P_OE", a_short_distances$int), ])
  PP=rbind(PP,d_PP)
  POE=rbind(POE,d_POE)
}

N_short_interactions <- short_distances[short_distances$int=="Total",]
N_short_interactions <- N_short_interactions[,c("sample","breaks","value")]
shorting_names <- c("Nut.0h","Nut.0h_1","Nut.0h_2","Nut.1h","Nut.1h_1","Nut.1h_2","Nut.10h","Nut.10h_1","Nut.10h_2")
names(shorting_names) <- c("DMSO_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed","DMSO_WT_BR1_recalibrated2WTDMSO_cutoff_5.ibed","DMSO_WT_BR2_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin1h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin1h_WT_BR1_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin1h_WT_BR2_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin10h_WT_merged_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin10h_WT_BR1_recalibrated2WTDMSO_cutoff_5.ibed","Nutlin10h_WT_BR2_recalibrated2WTDMSO_cutoff_5.ibed")
colors <- c("#fde725ff","#fde725ff","#fde725ff","#37b578ff","#37b578ff","#37b578ff","#43377fff","#43377fff","#43377fff")
names(colors) <- c("Nut.0h","Nut.0h_1","Nut.0h_2","Nut.1h","Nut.1h_1","Nut.1h_2","Nut.10h","Nut.10h_1","Nut.10h_2")

N_short_interactions$short_names <- shorting_names[match(N_short_interactions$sample,names(shorting_names))]
N_short_interactions$colors <- colors[match(N_short_interactions$short_names,names(colors))]

N_short_interactions$short_names <- factor(N_short_interactions$short_names, levels = c("Nut.0h","Nut.0h_1","Nut.0h_2","Nut.1h","Nut.1h_1", "Nut.1h_2","Nut.10h","Nut.10h_1", "Nut.10h_2"))

BR_interactions <- data.frame(N_short_interactions[!N_short_interactions$sample %like% "merge", ])
merged_interactions <- data.frame(N_short_interactions[N_short_interactions$sample %like% "merge", ])


N_short_interactions.df <- dcast(data = N_short_interactions,formula = breaks~short_names,fun.aggregate = sum,value.var = "value")
rownames(N_short_interactions.df) <- N_short_interactions.df$breaks

N_short_interactions.df_BR <- N_short_interactions.df[ , grepl( "_1|_2" , names( N_short_interactions.df ) ) ]
N_short_interactions.df_merged <- N_short_interactions.df[ , !grepl( "_1|_2" , names( N_short_interactions.df ) ) ]

samples <- list(N_short_interactions.df_BR[,-1],N_short_interactions.df_merged[,-1])

chisq <- chisq.test(samples[[2]])
corrplot(chisq$residuals, is.cor = FALSE,tl.col="black", method="circle",col=colorRampPalette(c("#23447f","#ffffff","#881010"))(100), sig.level = 0.005, insig = "blank") 



## Figure 4F ----

pkm <- load_interactions(file = "../../data/ComplementaryData/PCHiC/peakmatrix_merged_recalibrated2WTDMSO.txt") # chicago output
pkm2 <- as.matrix(pkm@elementMetadata[,6:c(ncol(pkm@elementMetadata)-1)])
pkm2 <-  asinh(pkm2)

scaled <- T
if (scaled) 
{
  madn <- mad(pkm2)
  mediann <- median(pkm2)
  pkm2[pkm2>mediann+3*madn] <- mediann+3*madn
}
pkm2 <- as.data.frame(pkm2)

# binarise peakmatrix
binarised <- as.data.frame(ifelse(pkm2>=asinh(5),yes=1,no=0))
binarised <- binarised %>% unite(group,CS_DMSO,CS_Nutlin3a_1h,CS_Nutlin3a_10h)
pkm2$group <- binarised$group

order <- as.data.frame(table(pkm2$group)) %>% arrange(desc(Freq))
order$Var1 <- as.character(order$Var1)

pkm2 <- pkm2 %>% arrange(factor(group,levels = order$Var1))

downsample <- as.numeric(0.1)
message("Selecting interactions based on downsample")
for(cl in unique(pkm2$group)){
  print(cl)
  sbst <- pkm2[pkm2$group==cl, ]
  
  num_sbst <- round(nrow(sbst)*downsample)
  set.seed(103) 
  
  if(nrow(sbst)>num_sbst){
    sbst <- sbst[sample(1:nrow(sbst), num_sbst), ]
  }
  
  if(cl==unique(pkm2$group)[1]){
    
    d = dist(sbst[,-ncol(sbst)])
    h = hclust(d, method="ward.D2")
    h$order
    newtbl <- sbst[h$order,]
    
  }else{ 
    d = dist(sbst[,-ncol(sbst)])
    h = hclust(d, method="ward.D2")
    h$order
    newtbl <- rbind(newtbl,sbst[h$order,]) 
  }
  
}

newtbl <- data.frame(newtbl,row.names = NULL)

col <- circlize::colorRamp2(breaks = c(0,asinh(4.99), asinh(5), mediann+3*madn), colors = c("#23447f","#91A1BF", "#C38787", "#881010"))


# visualize
x <- ComplexHeatmap::Heatmap(as.matrix(Fig_4_results$newtbl[,-c(4)]),col=col,cluster_rows = T,cluster_columns = F,split = Fig_4_results$newtbl$group,border=T,show_heatmap_legend = F)
ComplexHeatmap::draw(x)

pdf(paste0("heatmap.pdf"),width = 9,height = 12)
grid.newpage()
ComplexHeatmap::draw(ComplexHeatmap::Legend(col_fun=col))
dev.off()
