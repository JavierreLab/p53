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

pca_res <- prcomp(t(pkm2))
res <- as.data.frame(pca_res$x)
res$sample <- rownames(res)
res$condition <- gsub("_BR[1-2]","",res$sample)

# visualize
ggplot(res,aes(PC1,PC2,col=condition))+ 
  geom_point(size=5)+
  theme_classic2() +
  scale_color_manual(values=c("#fde725ff","#43377fff","#37b578ff"),label=c("Nut 0h","Nut 10h","Nut 1h")) +
  scale_fill_manual(values=c("#fde725ff","#43377fff","#37b578ff")) +
  xlab(paste0("PC1 (",summary(pca_res)$importance[2,1]*100,"%)")) +
  ylab(paste0("PC2 (",summary(pca_res)$importance[2,2]*100,"%)")) + 
  labs(fill="",color="") +
  theme(legend.position = "top")


## Figure 4B ----

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


## Figure 4C ----

### different script


## Figure 4D ----

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
x <- ComplexHeatmap::Heatmap(as.matrix(newtbl[,-c(4)]),col=col,cluster_rows = T,cluster_columns = F,split = newtbl$group,border=T,show_heatmap_legend = F)
ComplexHeatmap::draw(x)

pdf(paste0("heatmap.pdf"),width = 9,height = 12)
grid.newpage()
ComplexHeatmap::draw(ComplexHeatmap::Legend(col_fun=col))
dev.off()
