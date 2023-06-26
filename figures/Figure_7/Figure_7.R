library(tidyverse)




# Figure 7A ----

dds <- read_rds("../../data/ComplementaryData/H3K27ac/degron/normalised_peakmatrix_H3K27ac_WT_KD.Rds") # file too large for github
dds_rlog <- rlog(dds,blind=T)
dds_rlog <- assay(dds_rlog)

col <- circlize::colorRamp2(breaks = c(0.8,0.9,0.95,1),colors = c("#ffffff","#dee0eb","#bec1d8","#3b528b"))

# visualise
ComplexHeatmap::Heatmap(cor(dds_rlog),col=col,heatmap_legend_param = list(title = "Pearson Cor."))


# Figure 7B ----
pca_res <- prcomp(t(dds_rlog))
res <- as.data.frame(pca_res$x)
res$sample <- rownames(res)
res$condition <- gsub("_[1-3]","",res$sample)
res <- res %>% separate(condition,c("condition","timepoint"))
res$timepoint <- factor(res$timepoint, levels = c("Nutlin1h","DMSO"))

# visualise
ggplot(res,aes(PC1,PC2,col=timepoint,shape=condition))+  
  scale_color_manual(values=c("#37b578ff","#fde725ff"),guide=NULL) +
  scale_fill_manual(values=c("#37b578ff","#fde725ff"),guide=NULL) +
  scale_shape_manual(values=c(16,17),guide=NULL)+
  geom_point(size=5)+ 
  # ggforce::geom_mark_ellipse(aes(fill=condition)) +
  ggpubr::theme_classic2() + 
  xlab(paste0("PC1 (",summary(pca_res)$importance[2,1]*100,"%)")) + 
  ylab(paste0("PC2 (",summary(pca_res)$importance[2,2]*100,"%)"))


# Figure 7D ----

x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression.Rds")

x <- x %>% group_by(Cond,variable) %>% summarize(mean=mean(value),sd=sd(value))
x$Cond <- factor(x$Cond,levels = c("WTNut0h","KDNut0h","WTNut1h","KDNut1h","WTNut10h","KDNut10h"))

x <- x[!x$variable %in% c("PPP1CB","CHD4","ACTR8","TP53","TGFA","BAX","CDKN1A"),]
p<- ggplot(x, aes(x=variable, y=mean, fill=Cond)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  ggpubr::theme_classic2() + 
  scale_fill_manual(values = c("yellow","lightgrey","green","grey","purple","darkgrey"))

# visualise
p
