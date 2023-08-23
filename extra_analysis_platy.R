library(dartR)
library(fields)
library(ggplot2)
library(patchwork)
library(ape)
library(stringr)
library(dplyr)
library(tictoc)
library(data.table)

region_size <- 500000

file_path <- "/Users/s441489/final_platy/outflank"
list_files <- list.files(file_path,pattern = "26gen",full.names = TRUE)

list_files_tenter <- list.files(file_path,pattern = "tenter",full.names = TRUE)

file_path_pca <-  "/Users/s441489/final_platy/pca"
list_files_pca <- list.files(file_path_pca,full.names = TRUE,pattern = ".rda")

list_genes <- as.data.frame(matrix(nrow = 1,ncol =10 ))
colnames(list_genes) <- c("seqid","source","type","start","end","score","strand","phase","attributes", "length")

n_snps <- 0

scenario_fin <- as.data.frame(matrix(nrow = 1,ncol = 7))
colnames(scenario_fin) <- c( "position" ,"FST_tenter_coding","FST_tenter_noncoding", "FST_severn_coding","FST_severn_noncoding","scenario","chrom")

for(i in 1:length(list_files_pca)){
  
  chrom <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list_files_pca[i]))
  
  chrom <-  gsub(".*[_]", "\\1", chrom) 
  
  ########################GFF########################
  
  chrom_n <- as.numeric(gsub("([0-9]+).*$", "\\1", chrom))
  chrom_dir <- paste0("/Users/s441489/final_platy/chr_", chrom_n, "/")
  
  gff <- read.gff(paste0(chrom_dir,"chr_", chrom_n, "_mOrnAna1.pri.v4.gff"))
  
  genes <- gff[gff$type == "gene", ]
  
  genes_meta_tmp <- genes$attributes
  genes_meta <-
    as.data.frame(str_split(genes_meta_tmp, pattern = "=|;|-", simplify = TRUE))
  genes_meta <- genes_meta[, c(3, 5)]
  colnames(genes_meta) <- c("name", "ID")
  genes$attributes <- genes_meta$name
  genes$length <- genes$end - genes$start
  
  y_genes  <- list(genes$start, genes$end)
  
  ########################OUTFLANK########################
  t1 <- readRDS(list_files[i])
  t1_outlier <- t1$outflank$results
  t1_order <- t1_outlier
  t1_order$LocusName <- gsub(".*[_]([^_]+)[_].*", "\\1", t1_order$LocusName)
  t1_order$LocusName <- as.numeric(t1_order$LocusName)
  t1_order <- t1_order[order(t1_order$LocusName),]
  t1_order <- t1_order[!duplicated(t1_order$LocusName),]
  t1_order <- t1_order[complete.cases(t1_order),]
  t1_order <-   t1_order[,c("LocusName","He","FST","pvaluesRightTail")]
  t1_order$gene <- t1_order$LocusName %inrange% y_genes
  t1_order[t1_order$gene==FALSE,"gene"] <- "Non-coding Regions_Severn" 
  t1_order[t1_order$gene==TRUE,"gene"] <- "Coding Regions_Severn"
  t1_order$river <- "SEVERN"
  
  # TENTERFIELD 
  t1_t <- readRDS(list_files_tenter[i])
  t1_t_outlier <- t1_t$outflank$results
  t1_t_order <- t1_t_outlier
  t1_t_order$LocusName <- gsub(".*[_]([^_]+)[_].*", "\\1", t1_t_order$LocusName)
  t1_t_order$LocusName <- as.numeric(t1_t_order$LocusName)
  t1_t_order <- t1_t_order[order(t1_t_order$LocusName),]
  t1_t_order <- t1_t_order[!duplicated(t1_t_order$LocusName),]
  t1_t_order <- t1_t_order[complete.cases(t1_t_order),]
  t1_t_order <-   t1_t_order[,c("LocusName","He","FST","pvaluesRightTail")]
  t1_t_order$gene <- t1_t_order$LocusName %inrange% y_genes
  t1_t_order[t1_t_order$gene==FALSE,"gene"] <- "Non-coding Regions_Tenterfield"
  t1_t_order[t1_t_order$gene==TRUE,"gene"] <- "Coding Regions_Tenterfield"
  t1_t_order$river <- "TENTERFIELD"
  # 
  ########################PCA########################
  t1_pca <- readRDS(list_files_pca[i])
  t1_pca <- as.data.frame(t1_pca$loadings)
  t1_pca_order <- t1_pca
  t1_pca_order$LocusName <- rownames(t1_pca_order)
  t1_pca_order$LocusName <- gsub(".*[_]([^_]+)[_].*", "\\1",
                                 t1_pca_order$LocusName)
  t1_pca_order$LocusName <- as.numeric(t1_pca_order$LocusName)
  t1_pca_order <- t1_pca_order[order(t1_pca_order$LocusName),]
  t1_pca_order <- t1_pca_order[,c("LocusName","Axis1","Axis2","Axis3")]
  
  colnames(t1_pca_order) <- c(c("LocusName","PCA_Axis1","PCA_Axis2","PCA_Axis3"))
  
  join_df <- rbind(t1_order,t1_t_order)
  
  # join_df <- merge(join_df,t1_t_order,by="LocusName")
  
  join_df <- merge(join_df,t1_pca_order,by="LocusName")
  
  
  join_df_melt <- melt(join_df,id.vars=c("LocusName","gene","river"))
  
  pca_no_1 <- which(join_df_melt$variable == "PCA_Axis1" & join_df_melt$river == "TENTERFIELD")
  
  pca_no_2 <- which(join_df_melt$variable == "PCA_Axis2" & join_df_melt$river == "TENTERFIELD")
  
  pca_no_3 <- which(join_df_melt$variable == "PCA_Axis3" & join_df_melt$river == "TENTERFIELD")
  
  pca_no <- c(pca_no_1,pca_no_2,pca_no_3)
  
  join_df_melt <- join_df_melt[-pca_no,]
  
  n_snps <- n_snps + nrow(join_df)
  
  t1 <- join_df_melt[which(join_df_melt$variable=="FST"),]
  
  cols <- c("Coding Regions_Severn"= "Coding Regions_Severn","Non-coding Regions_Severn"="Non-coding Regions_Severn","Coding Regions_Tenterfield"="Coding Regions_Tenterfield","Non-coding Regions_Tenterfield"="Non-coding Regions_Tenterfield")
  
  p <- 
    ggplot(t1,aes(x=LocusName,y=value,color=gene))+
    stat_smooth(span = 1/100,n=round(tail(join_df_melt,1)$LocusName/region_size,0),se=FALSE,method = "loess") +
    scale_color_manual(values= cols)
  
  p2 <- ggplot_build(p)$data[[1]]
  
  # pvalue <- stats.bin(t1_order$LocusName,t1_order$pvaluesRightTail,breaks = seq(t1_order[1,"LocusName"],t1_order[nrow(t1_order),"LocusName"],10000))
  # # pvalue_sd <- sd(t1_order$pvaluesRightTail)
  # # pvalue_mean <- mean(t1_order$pvaluesRightTail)
  # pvalue_binned <- cbind(pvalue$stats[2,],pvalue$centers)
  # pvalue_binned <- as.data.frame(pvalue_binned)
  # colnames(pvalue_binned) <- c("pvalue","Position")
  # 
  # low_q <- quantile(pvalue_binned$pvalue,na.rm = TRUE, probs = seq(0, 1, 1 / 100),type=1)[2]
  # low_q_2 <-  pvalue_binned[which(pvalue_binned$pvalue<0.15),]
  # 
  # low_q_2$name <- NA
  # y  <- list(genes$start,genes$end)
  # for(x in 1:nrow(low_q_2)){
  #   genes_tmp <- which(low_q_2[x,"Position"] %between% y)
  #   list_genes <- rbind(list_genes,genes[genes_tmp,])
  #   if(length(genes_tmp)>0){
  #     low_q_2[x,"name"] <- genes[genes_tmp[1],"attributes"]
  #   }
  # }
  # 
  # low_q_2 <- low_q_2[complete.cases(low_q_2$name),]
  
  tenter_coding <- t1_t_order[which(t1_t_order$gene == "Coding Regions_Tenterfield" ),]
  tenter_noncoding <- t1_t_order[which(t1_t_order$gene == "Non-coding Regions_Tenterfield" ),]
  
  severn_coding <- t1_order[which(t1_order$gene == "Coding Regions_Severn" ),]
  severn_noncoding <- t1_order[which(t1_order$gene == "Non-coding Regions_Severn" ),]
  
  fst_tenter_coding_sd <- sd(tenter_coding$FST)
  fst_tenter_noncoding_sd <- sd(tenter_noncoding$FST)
  
  fst_severn_coding_sd <- sd(severn_coding$FST)
  fst_severn_noncoding_sd <- sd(severn_noncoding$FST)
  
  fst_tenter_coding_bin <- stats.bin(tenter_coding$LocusName,tenter_coding$FST,breaks = seq(tenter_coding[1,"LocusName"],tenter_coding[nrow(tenter_coding),"LocusName"],region_size))
  fst_tenter_coding_binned <- cbind(fst_tenter_coding_bin$stats[2,],fst_tenter_coding_bin$centers)
  fst_tenter_coding_binned <- as.data.frame(fst_tenter_coding_binned)
  colnames(fst_tenter_coding_binned) <- c("FST_tenter_coding","Position")
  
  fst_tenter_noncoding_bin <- stats.bin(tenter_noncoding$LocusName,tenter_noncoding$FST,breaks = seq(tenter_noncoding[1,"LocusName"],tenter_noncoding[nrow(tenter_noncoding),"LocusName"],region_size))
  fst_tenter_noncoding_binned <- cbind(fst_tenter_noncoding_bin$stats[2,],fst_tenter_noncoding_bin$centers)
  fst_tenter_noncoding_binned <- as.data.frame(fst_tenter_noncoding_binned)
  colnames(fst_tenter_noncoding_binned) <- c("FST_tenter_noncoding","Position")
  
  fst_severn_coding_bin <- stats.bin(severn_coding$LocusName,severn_coding$FST,breaks = seq(severn_coding[1,"LocusName"],severn_coding[nrow(severn_coding),"LocusName"],region_size))
  fst_severn_coding_binned <- cbind(fst_severn_coding_bin$stats[2,],fst_severn_coding_bin$centers)
  fst_severn_coding_binned <- as.data.frame(fst_severn_coding_binned)
  colnames(fst_severn_coding_binned) <- c("FST_severn_coding","Position")
  
  fst_severn_noncoding_bin <- stats.bin(severn_noncoding$LocusName,severn_noncoding$FST,breaks = seq(severn_noncoding[1,"LocusName"],severn_noncoding[nrow(severn_noncoding),"LocusName"],region_size))
  fst_severn_noncoding_binned <- cbind(fst_severn_noncoding_bin$stats[2,],fst_severn_noncoding_bin$centers)
  fst_severn_noncoding_binned <- as.data.frame(fst_severn_noncoding_binned)
  colnames(fst_severn_noncoding_binned) <- c("FST_severn_noncoding","Position")
  
  fst_merge <- cbind(fst_tenter_coding_binned$Position, 
                     fst_tenter_coding_binned$FST_tenter_coding,
                     fst_tenter_noncoding_binned$FST_tenter_noncoding,
                     fst_severn_coding_binned$FST_severn_coding,
                     fst_severn_noncoding_binned$FST_severn_noncoding
  )
  
  fst_merge <- as.data.frame(fst_merge)
  colnames(fst_merge) <- c("position",
                           "FST_tenter_coding",
                           "FST_tenter_noncoding",
                           "FST_severn_coding",
                           "FST_severn_noncoding"
  )
  fst_merge <- fst_merge[complete.cases(fst_merge),]
  
  #   a) Adaptive loci diverged, neutral loci did not diverged; 
  # - FST coding Severn (mean region 100Kbp) > 1 SD + FST coding Tenterfield  (mean region 100Kbp)
  # - FST coding Severn (mean region 100Kbp) > 1 SD + FST non-coding Severn  (mean region 100Kbp)
  
  scenario_a <- fst_merge[which(fst_merge$FST_severn_coding > 
                                  (fst_tenter_coding_sd + fst_merge$FST_tenter_coding) &
                                  fst_merge$FST_severn_coding > 
                                  (fst_severn_noncoding_sd + fst_merge$FST_severn_noncoding)),]
  if(nrow(scenario_a)>0){
    scenario_a$scenario <- "a"
    scenario_a$chrom <- chrom_n
  }else{
    scenario_a[1,] <- NA
    scenario_a$scenario <- "a"
    scenario_a$chrom <- chrom_n
  }
  
  # b) Neutral loci diverged, adaptive loci did not diverged; 
  # - FST non-coding Severn (mean region 100Kbp) > 1 SD + FST non-coding Tenterfield (mean region 100Kbp)
  # - FST non-coding Severn (mean region 100Kbp) > 1 SD + FST coding Severn (mean region 100Kbp)
  
  scenario_b <- fst_merge[which(fst_merge$FST_severn_noncoding > 
                                  (fst_tenter_noncoding_sd + fst_merge$FST_tenter_noncoding) &
                                  fst_merge$FST_severn_noncoding > 
                                  (fst_severn_coding_sd + fst_merge$FST_severn_coding)),]
  
  if(nrow(scenario_b)>0){
    scenario_b$scenario <- "b"
    scenario_b$chrom <- chrom_n
  }else{
    scenario_b[1,] <- NA
    scenario_b$scenario <- "b"
    scenario_b$chrom <- chrom_n
  }
  
  # c) Adaptive loci diverged, neutral loci diverged. 
  # - FST coding Severn (mean region 100Kbp) > 1 SD + FST coding Tenterfield (mean region 100Kbp)
  # - FST non-coding Severn (mean region 100Kbp) > 1 SD + FST non-coding Tenterfield  (mean region 100Kbp)
  
  scenario_c <- fst_merge[which(fst_merge$FST_severn_coding > 
                                  (fst_tenter_coding_sd + fst_merge$FST_tenter_coding) &
                                  fst_merge$FST_severn_noncoding > 
                                  (fst_tenter_noncoding_sd + fst_merge$FST_tenter_noncoding)),]
  
  if(nrow(scenario_c)>0){
    scenario_c$scenario <- "c"
    scenario_c$chrom <- chrom_n
  }else{
    scenario_c[1,] <- NA
    scenario_c$scenario <- "c"
    scenario_c$chrom <- chrom_n
  }
  
  scenario_tmp <- as.data.frame(rbind(scenario_a,scenario_b,scenario_c))
  
  scenario_tmp <- scenario_tmp[complete.cases(scenario_tmp),]
  
  scenario_fin <- rbind(scenario_fin,scenario_tmp)
  
  if(nrow(scenario_tmp)>0){
    
    join_df_melt$stat <- factor(join_df_melt$variable, levels=c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3"))
    
    # cols <- c("Coding Regions"= "dodgerblue","Non-coding Regions"="dodgerblue4","Coding Regions_Tenter"="gold","Non-coding Regions_Tenter"="gold4")
    
    cols <- c("Coding Regions_Severn"= "lightgreen","Non-coding Regions_Severn"="darkgreen","Coding Regions_Tenterfield"="orange","Non-coding Regions_Tenterfield"="orange4")
    # 
    # cols <- c("Coding Regions"= "deeppink","Non-coding Regions"="deepskyblue3")
    
    f_labels <- data.frame(stat = factor("FST",levels = c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3")),
                           x = scenario_tmp$position ,
                           label =  scenario_tmp$scenario)
    
    p1 <- 
      # ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene,fill=gene))+
      ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene))+
      facet_wrap(~stat,scales="free_y",ncol=1)+
      stat_smooth(span = 1/100,se=FALSE,method = "loess",n=round(tail(join_df_melt,1)$LocusName/region_size,0))+
      theme_classic() +
      scale_x_continuous(breaks=seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000),labels=as.character(round(seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000)/1000000))) +
      xlab("Genome position (Mbp)")+
      labs(title = paste("Chromosome",chrom,nrow(join_df),"SNPs"), color = "") +
      scale_color_manual(values= cols) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      geom_text(data=f_labels,aes(label = label),x =f_labels$x, y = -0.1,inherit.aes = FALSE,size =6) 
    
    ggsave(paste0("v3_platypus_26_chr_",chrom,".pdf"),  width = 10, height =10, units = "in", dpi="retina", bg = "transparent" )
    
  }else{
    
    join_df_melt$stat <- factor(join_df_melt$variable, levels=c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3"))
    
    # cols <- c("Coding Regions"= "dodgerblue","Non-coding Regions"="dodgerblue4","Coding Regions_Tenter"="gold","Non-coding Regions_Tenter"="gold4")
    
    cols <- c("Coding Regions_Severn"= "lightgreen","Non-coding Regions_Severn"="darkgreen","Coding Regions_Tenterfield"="orange","Non-coding Regions_Tenterfield"="orange4")
    # 
    # cols <- c("Coding Regions"= "deeppink","Non-coding Regions"="deepskyblue3")
    
    
    p1 <- 
      # ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene,fill=gene))+
      ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene))+
      facet_wrap(~stat,scales="free_y",ncol=1)+
      stat_smooth(span = 1/100,se=FALSE,method = "loess")+
      theme_classic() +
      scale_x_continuous(breaks=seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000),labels=as.character(round(seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000)/1000000))) +
      xlab("Genome position (Mbp)")+
      labs(title = paste("Chromosome",chrom,nrow(join_df),"SNPs"), color = "") +
      scale_color_manual(values= cols) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5)) 
    # +
    #   geom_text(data=f_labels,aes(label = label),x =f_labels$x, y = -0.1,inherit.aes = FALSE,size =6) 
    
    p2 <- ggplot_build(p1)
    
    head(p2$data[[2]])
    
    ggsave(paste0("v3_platypus_26_chr_",chrom,".pdf"),  width = 10, height =10, units = "in", dpi="retina", bg = "transparent" )
    
  }
  
}
write.csv(scenario_fin,file = "scenario_fin.csv")
print(n_snps)

# scenario_a = 104
# scenario_b = 44
# scenario_a = 193



# write.csv(list_genes,file = "list_genes.csv")
