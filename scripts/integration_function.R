
p53_genes_interactions <- function(ints, p53.chip, promoter.annot, at.promoter, interaction.type){
  
  if(at.promoter){
    ints <- interactionsByRegions(ints,p53.chip)
    df <- data.frame(bait=as.character(anchorOne(ints)),
                     oe=as.character(anchorTwo(ints)),
                     gene_I=ints$gene_I,
                     gene_II=ints$gene_II,
                     int=ints$int,
                     overlap_I=ints$overlap_I,
                     overlap_II=ints$overlap_II)
    
    df$overlapping_region <- ifelse(df$overlap_I,yes=df$bait,no=df$oe)
    df$other_region <- ifelse(df$overlap_I,yes=df$oe,no=df$bait)
    
    overlapping_regions <- makeGRangesFromDataFrame(df %>% separate(overlapping_region,c("chr","start","end"),sep=":|-"),keep.extra.columns = T)
    merge <- mergeByOverlaps(p53.chip,overlapping_regions)
    merge <- data.frame(p53=as.character(merge$p53.chip),overlapping_region=as.character(merge$overlapping_regions))
    df <- left_join(df,merge)
    
    df$interacting_gene <- ifelse(df$overlap_I,yes=df$gene_II,no=df$gene_I)
    df <- df %>% separate_rows(interacting_gene,sep = ",")
    df$gene_I <- NULL
    df$gene_II <- NULL
    
    merge <- mergeByOverlaps(p53.chip,promoter.annot)
    aux <- data.frame(p53=as.character(merge$p53.chip),target_gene=merge$ensembl_gene_id)
    df <- unique(left_join(df,aux))
    return(df)
  }else if(!at.promoter){
    if(interaction.type=="P_OE"){
      
      # integrate PCHiC with p53 
      ints <- HiCaptuRe::interactionsByRegions(ints,p53.chip)
      
      
      # keep subset PCHiC info as dataframe
      df <- data.frame(gene_I=ints@elementMetadata$gene_I,
                       gene_II=ints@elementMetadata$gene_II,
                       overlap_I=ints@elementMetadata$overlap_I,
                       overlap_II=ints@elementMetadata$overlap_II,
                       bait=as.character(anchorOne(ints)),
                       oe=as.character(anchorTwo(ints)),
                       int_type=ints$int)
      df <- df %>% group_by(int_type) %>% group_split()
      
      df2 <- lapply(df,function(x){
        if(unique(x$int_type) %in% c("OE_P_I","P_OE_II")){
          x <- separate_rows(x,gene_I,sep = ",")
          x <- x %>% separate(oe, c("chr","start","end"),sep = ":|-",remove=F)
          oe <- makeGRangesFromDataFrame(x,keep.extra.columns = T)
          merge <- mergeByOverlaps(oe,p53.chip)
          link <- unique(data.frame(p53=as.character(merge$p53.chip),ensg=merge$gene_I,bait=merge$bait,oe=as.character(merge$oe),overlap_I=merge$overlap_I,overlap_II=merge$overlap_II,int_type=merge$int_type))
          return(link)
        }else if(unique(x$int_type) %in% c("P_P_I","P_P_II")){
          x <- separate_rows(x,gene_I,sep = ",")
          x <- separate_rows(x,gene_II,sep = ",")
          
          x$ensg <- ifelse(x$overlap_I & !x$overlap_II,yes=x$gene_II,no=ifelse(!x$overlap_I & x$overlap_II,yes=x$gene_I,no=paste0(x$gene_I,",",x$gene_II)))
          x <- separate_rows(x,ensg,sep = ",")
          x <- unique(x[,!colnames(x) %in% c("gene_I","gene_II")])
          
          
          x$regions <- ifelse(x$overlap_I & !x$overlap_II,yes=x$bait,no=ifelse(!x$overlap_I & x$overlap_II,yes=x$oe,no=paste0(x$bait,",",x$oe)))
          x <- separate_rows(x,regions,sep = ",")
          
          x <- x %>% separate(regions, c("chr","start","end"),sep = ":|-",remove=F)
          regions <- makeGRangesFromDataFrame(x,keep.extra.columns = T)
          merge <- mergeByOverlaps(regions,p53.chip)
          merge <- merge[merge$p53.chip %in% subsetByOverlaps(merge$p53.chip,promoter.annot,invert = T),]
          link <- unique(data.frame(p53=as.character(merge$p53.chip),ensg=merge$ensg,bait=merge$bait,oe=as.character(merge$oe),overlap_I=merge$overlap_I,overlap_II=merge$overlap_II,int_type=merge$int_type))
        }
      })
      df2 <- do.call("rbind",df2)
      return(df2)
      
      
    }else if(interaction.type=="P_P"){

    }
    
  }
  
  return(df)
}

