library(GenomicInteractions)
library(HiCaptuRe)
library(tidyverse)
library(DESeq2)

setwd("~/MEGAsync/scripts/p53/2_SUBMISSION/github/figures/Figure_7/")

# Figure 7A ----
pkm <- HiCaptuRe::load_interactions("../../data/ComplementaryData/PCHiC/peakmatrix_merged_Nut1h_WT_KD_recalibrated2WTDMSO.txt")
aux <- data.frame(bait=as.character(anchorOne(pkm)),
                  oe=as.character(anchorTwo(pkm)),
                  chiscore_WT=(pkm@elementMetadata$CS_Nutlin3a_1h),
                  chiscore_KD=(pkm@elementMetadata$CS_Rad21KD_Nutlin3a_1h))

df <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1h.Rds")
df <- unique(df[,c(1:13)])
df <- left_join(df,aux)
df <- unique(df[,c("chiscore_WT","chiscore_KD")]) 
df <- reshape2::melt(df)
df$variable <- as.character(df$variable)
df$variable <- factor(ifelse(df$variable %in% "chiscore_WT",yes="WT",no="KD"), levels = c("WT","KD"))

# visualise
ggplot(to.plot,aes(variable,value,fill=variable)) + 
  geom_violin(aes(col=variable),fill=NA)+
  geom_boxplot(width=0.3,outlier.shape = NA) + 
  ggpubr::theme_classic2() +
  scale_fill_manual(values=c("#37b578ff","darkgrey"),guide=NULL) +
  scale_color_manual(values=c("#37b578ff","darkgrey"),guide=NULL) +
  labs(x="Condition at Nut1h",y="Chicago Score of p53 - distal target gene interactions",fill="") +
  ylim(0,15)


# Figure 7B ----

df <- read_rds("../../data/ComplementaryData/PCHiC_integration/distal_integration_Nutlin1h.Rds")

ints_kd <- HiCaptuRe::load_interactions("../../data/ComplementaryData/PCHiC/Nutlin1h_KD_merged_recalibrated2WTDMSO_cutoff_5.ibed")
ints_kd <- data.frame(bait=as.character(anchorOne(ints_kd)),
                     oe=as.character(anchorTwo(ints_kd)),
                     chiscore=ints_kd@elementMetadata$CS)

df$interaction_KD <- paste0(df$bait,",",df$oe) %in% unique(paste0(ints_kd$bait,",",ints_kd$oe))

to.plot <- unique(df[,c("bait","oe","interaction_KD")]) %>% group_by(interaction_KD) %>% summarise(KD=length(interaction_KD))
to.plot$WT <- c(0,sum(to.plot$KD))
to.plot <- reshape2::melt(to.plot)                                                                                        
to.plot$interaction_KD <- ifelse(to.plot$interaction_KD,"Maint.","Lost")
to.plot <- to.plot[to.plot$variable=="KD",]

# visualise
ggplot(to.plot,aes(variable,value,fill=interaction_KD)) + 
  geom_col(width=0.7,col="black") + 
  ggpubr::theme_classic2() +
  scale_fill_manual(values=c("lightgrey","darkgrey")) +
  labs(x="Interactions \n at Nut 1h",y="Number of interactions",fill="",title="") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "top")


# Figure 7C ----
to.plot <- unique(df[,c("genomic_distance","interaction_KD")]) 
to.plot$interaction_KD <- as.logical(to.plot$interaction_KD)
to.plot$interaction_KD <- ifelse(to.plot$interaction_KD,"Maint.","Lost")
to.plot <- to.plot %>% group_by(interaction_KD) %>% summarize(med_distance=median(distance_p53_ensg,na.rm=T),distance_p53_ensg=distance_p53_ensg)

# visualise
ggplot(to.plot,aes(interaction_KD,log(genomic_distance),fill=interaction_KD)) + 
  geom_violin(aes(col=interaction_KD),width=0.5,fill=NA) +
  geom_boxplot(width=0.2,outlier.shape = NA) + 
  ggpubr::theme_classic2() +
  scale_colour_manual(values=c("lightgrey","darkgrey"),guide=NULL) +
  scale_fill_manual(values=c("lightgrey","darkgrey"),guide=NULL) +
  labs(x="Interaction dynamics",y="Genomic distance log(bp)",fill="")


# Figure 7E ----
x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression_distal_genes_WT_KD.Rds")
x <- x[,-c(11:12)]
x <- reshape2::melt(x,colnames(x)[1:4])
x$variable <- gsub("RelativeExpression_BR","",as.character(x$variable))
colnames(x)[5] <- "BR"
x <- reshape2::dcast(x,formula = `Cell.type/Cell.line`+Condition+Timepoint+BR~GeneID,value.var = "value")
x <- x[,match(c("Cell.type/Cell.line","Condition","Timepoint","BR","TGFBR2","GSR","LAMC1","PLK2","JAG2","BRD7","PPM1D","TENT4B","TEP1","CD9","S100A1","PCM1"),colnames(x))]

step1 <- x[x$Condition=="WT",]
y0h <- c("WT","Nut0h",colMeans(step1[step1$Timepoint=="Nut0h",5:ncol(step1)],na.rm = T))
y1h <- c("WT","Nut1h",colMeans(step1[step1$Timepoint=="Nut1h",5:ncol(step1)],na.rm = T))
y10h <- c("WT","Nut10h",colMeans(step1[step1$Timepoint=="Nut10h",5:ncol(step1)],na.rm = T))
y <- as.data.frame(rbind(y0h,y1h,y10h))
colnames(y)[1:2] <- c("Condition","Timepoint")
rownames(y) <- 1:3
step1 <- y

step2 <- x
for(i in 5:ncol(step2)){
  gene <- colnames(step2)[i]
  for(j in 1:nrow(step2)){
    aux_timepoint <- step2[j,]$Timepoint
    step2[j,i] <- step2[j,i]/as.numeric(step1[step1$Timepoint==aux_timepoint,colnames(step1)==gene])
  }
  
}
x <- step2

to.plot <- reshape2::melt(x,colnames(x)[1:4])
to.plot <- to.plot[!is.na(to.plot$value),]
to.plot$variable <- NULL

to.plot$Timepoint_Condition <- paste0(to.plot$Condition,to.plot$Timepoint)
to.plot$Timepoint <- factor(to.plot$Timepoint,levels = c(unique(grep("Nut0h",to.plot$Timepoint,value = T)),
                                                         unique(grep("Nut1h",to.plot$Timepoint,value = T)),
                                                         unique(grep("Nut10h",to.plot$Timepoint,value = T))))

to.plot$Condition <- factor(to.plot$Condition,levels=c("WT","KD"))
to.plot$Timepoint_Condition <- factor(to.plot$Timepoint_Condition,levels=c("WTNut0h","KDNut0h",
                                                                           "WTNut1h","KDNut1h",
                                                                           "WTNut10h","KDNut10h"))

p<- ggplot(to.plot, aes(x=Condition, y=value,fill=Timepoint_Condition)) +
  facet_grid(~Timepoint) +
  geom_violin(aes(col=Timepoint_Condition),width=0.9,fill=NA,position=position_dodge(width = 0.9)) +
  geom_boxplot(width=0.6,outlier.shape = NA,col="black",position=position_dodge(width = 0.9)) +
  ggpubr::theme_classic2() +
  scale_fill_manual(values = c("#fce8267f","grey","#37b578ff","grey","#43377fff","grey"),guide=NULL) +
  scale_color_manual(values = c("#fce8267f","grey","#37b578ff","grey","#43377fff","grey"),guide=NULL) +
  labs(x="Condition",y="WT-normalized nascent \n mRNA expression (FC)") 
p

# Statistics
x <- reshape2::melt(x,colnames(x)[1:4])
x <- x[!is.na(x$value),]
j <- 1
res <- list()
for(timepoint in unique(x$Timepoint)){
  aux.test <- wilcox.test(x[x$Timepoint==timepoint & x$Condition=="WT",]$value,x[x$Timepoint==timepoint & x$Condition=="KD",]$value)
  res[[j]] <- c(timepoint,aux.test[[3]])
  j <- j+1
}
res


# Figure 7F ----

x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression_distal_genes_WT_KD.Rds")
x <- x[c(1:4,11:12)]
x$Condition <- paste0(x$Condition,x$Timepoint)
x$Condition <- factor(x$Condition,levels = c("WTNut0h","KDNut0h","WTNut1h","KDNut1h","WTNut10h","KDNut10h"))
x$GeneID <- factor(x$GeneID,c("TGFBR2","GSR","LAMC1","PLK2","JAG2","BRD7","PPM1D","TENT4B","TEP1","CD9","S100A1","PCM1"))

p<- ggplot(x, aes(x=GeneID, y=Mean, fill=Condition)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-StandardDeviation,ymax=Mean+StandardDeviation), width=.2,
                position=position_dodge(.9)) + 
  ggpubr::theme_classic2() + 
  scale_fill_manual(values = c("yellow","lightgrey","green","grey","purple","darkgrey"))+
  ggbreak::scale_y_cut(breaks = c(8.5)) +
  labs(x="",t="Realtive nascent mRNA expression (FC)")

# visualise
p

# Statistics
x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression_distal_genes_WT_KD.Rds")
x <- reshape2::melt(x[,-c(11:12)],colnames(x)[1:4])
x <- x[!is.na(x$value),]


df_list <- list()
j <- 1
for(timepoint in grep("Nut0h",unique(x$Timepoint),value=T,invert=T)){
  df <- as.data.frame(matrix(nrow=4))[,-1]
  
  for(gene in unique(x[x$Timepoint==timepoint, ]$GeneID)){
    print(paste0(timepoint,": ",gene))
    res <- t.test(x=x[x$Timepoint==timepoint & x$GeneID==gene & x$Condition=="WT", ]$value,
                  y=x[x$Timepoint==timepoint & x$GeneID==gene & x$Condition=="KD", ]$value,
                  alternative = "greater")
    pvalue <- res[[3]]
    vec <- unlist(c(gene, timepoint, pvalue))
    df <- rbind(df,as.data.frame(t(as.data.frame(vec))))
  }
  rownames(df) <- 1:nrow(df)
  colnames(df) <- c("variable","Timepoint","pvalue")
  df_list[[j]] <- df
  j<-j+1
}
df_stats <- do.call("rbind",df_list)
df_stats$pvalue <- ifelse(as.numeric(df_stats$pvalue)<0.05,yes="*",no="ns")
df_stats$pvalue[is.na(df_stats$pvalue)] <- ""
df_stats$Condition <- "WT vs. KD"

# Figure 7I ----
x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression_CRISPR_validation.Rds")
x <- x[x$GeneID %in% c("S100A1","PCM1"),]
x <- x[,-c(5:10)]
x$Condition <- gsub("clone ","",x$Condition)

x1 <- x[x$GeneID=="S100A1",]
x2 <- x[x$GeneID=="PCM1",]

x1$Condition <- factor(x1$Condition,levels = c("WT","p53BS-/- S100A1 1","p53BS-/- S100A1 2"))
p1<- ggplot(x1, aes(x=Timepoint, y=Mean, fill=Condition)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-StandardDeviation,ymax=Mean+StandardDeviation), width=.2,
                position=position_dodge(.9)) + 
  ggpubr::theme_classic2() + 
  scale_fill_manual(values = c("white","darkgrey","black"))+
  labs(x="Timepoint",y="Realtive nascent mRNA expression (FC)")


x2$Condition <- factor(x2$Condition,levels = c("WT","p53BS-/- PCM1 1","p53BS-/- PCM1 2"))
p2<- ggplot(x2, aes(x=Timepoint, y=Mean, fill=Condition)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean-StandardDeviation,ymax=Mean+StandardDeviation), width=.2,
                position=position_dodge(.9)) + 
  ggpubr::theme_classic2() + 
  scale_fill_manual(values = c("white","darkgrey","black"))+
  labs(x="Timepoint",y="Realtive nascent mRNA expression (FC)")

# visualise
ggarrange(p1,p2)


# statistics
x <- read_rds("../../data/ComplementaryData/qrtpcr/qRTPCR_relative_expression_CRISPR_validation.Rds")
x <- x[,-c(11:12)]
x <- reshape2::melt(x,colnames(x)[1:4])
x <- x[!is.na(x$value),]
x <- x[x$GeneID %in% c("S100A1","PCM1"),]

df_list <- list()
j <- 1
for(condition in grep("WT",unique(x$Condition),invert=T,value=T)){
  df <- as.data.frame(matrix(nrow=4))[,-1]
  
  for(gene in unique(x[x$Condition==condition, ]$GeneID)){
    print(paste0(condition,": ",gene))
    res <- t.test(x=x[x$Condition=="WT" & x$Timepoint=="Nut1h" & x$GeneID==gene, ]$value,
                  y=x[x$Condition==condition & x$Timepoint=="Nut1h" &  x$GeneID==gene, ]$value,
                  alternative="greater")
    pvalue <- res[[3]]
    vec <- unlist(c(gene, condition,"Nut1h", pvalue))
    df <- rbind(df,as.data.frame(t(as.data.frame(vec))))
    
    print(paste0(condition,": ",gene))
    res <- t.test(x=x[x$Condition=="WT" & x$Timepoint=="Nut0h" & x$GeneID==gene, ]$value,
                  y=x[x$Condition==condition & x$Timepoint=="Nut0h" &  x$GeneID==gene, ]$value,
                  alternative="greater")
    pvalue <- res[[3]]
    vec <- unlist(c(gene, condition,"Nut0h", pvalue))
    df <- rbind(df,as.data.frame(t(as.data.frame(vec))))
    
  }
  rownames(df) <- 1:nrow(df)
  colnames(df) <- c("variable","Condition","Timepoint","pvalue")
  df_list[[j]] <- df
  j<-j+1
}

df_stats <- do.call("rbind",df_list)
df_stats$pvalue <- ifelse(as.numeric(df_stats$pvalue)<0.05,yes="*",no="ns")
df_stats$pvalue[is.na(df_stats$pvalue)] <- ""
df_stats$Condition <- paste0("WT vs. ",df_stats$Condition)
df_stats
