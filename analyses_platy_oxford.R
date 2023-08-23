library(dartR)
library(fields)
library(ggplot2)
library(patchwork)
library(ape)
library(stringr)
library(dplyr)
library(tictoc)
library(data.table)

region_size <- 100000
n_sd <- 2
span <- 200

file_path <- "/Users/s441489/final_platy/outflank"
list_files <- list.files(file_path,pattern = "ox_",full.names = TRUE)

# list_files_tenter <- list.files(file_path,pattern = "tenter",full.names = TRUE)

file_path_pca <-  "/Users/s441489/final_platy/pca"
list_files_pca <- list.files(file_path_pca,pattern = "ox_pca",full.names = TRUE)

list_genes <- as.data.frame(matrix(nrow = 1,ncol =10 ))
colnames(list_genes) <- c("seqid","source","type","start","end","score","strand","phase","attributes", "length")

n_snps <- 0

scenario_fin <- as.data.frame(matrix(nrow = 1,ncol = 7))
colnames(scenario_fin) <- c( "position" ,"FST_tenter_coding","FST_tenter_noncoding", "FST_severn_coding","FST_severn_noncoding","scenario","chrom")

gene_names <- as.data.frame(matrix(nrow =1,ncol = 39))
colnames(gene_names)<-c("seqid","source"   ,  "type", "start","end", 
                        "score","strand" ,    "phase","attributes", "length"    ,
                        "V1",   "V2",   "V3",   "V4",   "V5",  
                        "V6",   "V7",   "V8",   "V9",   "V10", 
                        "V11",  "V12",  "V13",  "V14",  "V15", 
                        "V16",  "V17",  "V18",  "V19",  "V20", 
                        "V21",  "V22",  "V23",  "V24",  "V25", 
                        "V26",  "V27",  "V28",  "V29"   )
gene_names$chr <- NA

chrom_platypus <- read.csv(paste0(getwd(), "/info/chrom_platypus.csv"))

# colnames(gene_names) <- c("seqid","source","type" ,"start","end", "score",  "strand", "phase","attributes", "length","chr" )

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
  
  mRNA <- gff[gff$type == "mRNA", ]
  mRNA_meta_tmp <- mRNA$attributes
  mRNA_meta <-
    as.data.frame(str_split(mRNA_meta_tmp, pattern = "=|;|-", simplify = TRUE))
  
  mRNA$length <- mRNA$end - mRNA$start
  mRNA <- cbind(mRNA,mRNA_meta)
  
  y_mRNA  <- list(mRNA$start, mRNA$end)
  
  coding <- genes[,c("start","end")]
  coding$region <- "Coding Regions"
  
  non_coding <- cbind(coding[,"end"]+1,coding[2:nrow(coding),"start"]-1)
  non_coding <- rbind(c(1,coding[1,"start"]-1),non_coding)
  non_coding[nrow(non_coding),2] <- gff[1,"end"]
  non_coding <- as.data.frame(non_coding)
  non_coding$region <- "Non-coding Regions"
  colnames(non_coding) <- c("start","end","region")
  non_coding$y <- 0
  
  tmp <- which(chrom_platypus$Name_b==chrom_n)
  
  last_snp_cod <- as.numeric(chrom_platypus[tmp,"end"])
  last_snp_cod <- which(non_coding$end < last_snp_cod)
  last_snp_cod <- tail(last_snp_cod,1)
  non_coding <- non_coding[1:last_snp_cod,]
  non_coding$length <- non_coding$end - non_coding$start
  #check whether is true that not real non-coding regions have negative value
  non_coding <- non_coding[which(non_coding$length>0),]
  
  last_snp_non <- as.numeric(chrom_platypus[tmp,"end"])
  last_snp_non <- which(coding$end < last_snp_non)
  last_snp_non <- tail(last_snp_non,1)
  coding <- coding[1:last_snp_non,]
  coding$y <- -Inf
  coding$length <- coding$end - coding$start
  
  regions <- rbind(coding,non_coding)
  
  ########################OUTFLANK########################
  t1 <- readRDS(list_files[i])
  t1_outlier <- t1$results
  t1_order <- t1_outlier
  t1_order$LocusName <- gsub(".*[_]([^_]+)[_].*", "\\1", t1_order$LocusName)
  t1_order$LocusName <- as.numeric(t1_order$LocusName)
  t1_order <- t1_order[order(t1_order$LocusName),]
  t1_order <- t1_order[!duplicated(t1_order$LocusName),]
  t1_order <- t1_order[complete.cases(t1_order),]
  t1_order <-   t1_order[,c("LocusName","He","FST","pvaluesRightTail")]
  t1_order$gene <- t1_order$LocusName %inrange% y_genes
  # t1_order[t1_order$gene==FALSE,"gene"] <- "Non-coding Regions_Severn" 
  # t1_order[t1_order$gene==TRUE,"gene"] <- "Coding Regions_Severn"
  t1_order[t1_order$gene==FALSE,"gene"] <- "Non-coding Regions" 
  t1_order[t1_order$gene==TRUE,"gene"] <- "Coding Regions"
  # t1_order$river <- "SEVERN"
  
  # TENTERFIELD 
  # t1_t <- readRDS(list_files_tenter[i])
  # t1_t_outlier <- t1_t$outflank$results
  # t1_t_order <- t1_t_outlier
  # t1_t_order$LocusName <- gsub(".*[_]([^_]+)[_].*", "\\1", t1_t_order$LocusName)
  # t1_t_order$LocusName <- as.numeric(t1_t_order$LocusName)
  # t1_t_order <- t1_t_order[order(t1_t_order$LocusName),]
  # t1_t_order <- t1_t_order[!duplicated(t1_t_order$LocusName),]
  # t1_t_order <- t1_t_order[complete.cases(t1_t_order),]
  # t1_t_order <-   t1_t_order[,c("LocusName","He","FST","pvaluesRightTail")]
  # t1_t_order$gene <- t1_t_order$LocusName %inrange% y_genes
  # t1_t_order[t1_t_order$gene==FALSE,"gene"] <- "Non-coding Regions_Tenterfield"
  # t1_t_order[t1_t_order$gene==TRUE,"gene"] <- "Coding Regions_Tenterfield"
  # t1_t_order$river <- "TENTERFIELD"
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
  
  # join_df <- rbind(t1_order,t1_t_order)
  
  # join_df <- merge(join_df,t1_t_order,by="LocusName")
  
  # join_df <- merge(join_df,t1_pca_order,by="LocusName")
  
  join_df <- merge(t1_order,t1_pca_order,by="LocusName")
  
  #join_df_melt <- melt(join_df,id.vars=c("LocusName","gene","river"))
  
  join_df_melt <- melt(join_df,id.vars=c("LocusName","gene"))
  
  # pca_no_1 <- which(join_df_melt$variable == "PCA_Axis1" & join_df_melt$river == "TENTERFIELD")
  # 
  # pca_no_2 <- which(join_df_melt$variable == "PCA_Axis2" & join_df_melt$river == "TENTERFIELD")
  # 
  # pca_no_3 <- which(join_df_melt$variable == "PCA_Axis3" & join_df_melt$river == "TENTERFIELD")
  # 
  # pca_no <- c(pca_no_1,pca_no_2,pca_no_3)
  # 
  # join_df_melt <- join_df_melt[-pca_no,]
  # 
  n_snps <- n_snps + nrow(join_df)
  
  t1 <- join_df_melt[which(join_df_melt$variable=="FST"),]
  
  # cols <- c("Coding Regions_Severn"= "Coding Regions_Severn","Non-coding Regions_Severn"="Non-coding Regions_Severn","Coding Regions_Tenterfield"="Coding Regions_Tenterfield","Non-coding Regions_Tenterfield"="Non-coding Regions_Tenterfield")
  
  cols <- c("Coding Regions"= "Coding Regions","Non-coding Regions"="Non-coding Regions")
  
  span_2 <- nrow(join_df)/2
  
  
  p <- 
    ggplot(t1,aes(x=LocusName,y=value,color=gene))+
    stat_smooth(span = span/ span_2,n=round(tail(join_df_melt,1)$LocusName/region_size,0),se=FALSE,method = "loess") +
    scale_color_manual(values= cols)
  
  p2 <- ggplot_build(p)$data[[1]]
  
coding <- p2[which(p2$colour == "Coding Regions" ),]
noncoding <- p2[which(p2$colour == "Non-coding Regions" ),]
  
coding_sd <- sd(coding$y)
noncoding_sd <- sd(noncoding$y)
  
  # fst_tenter_coding_bin <- stats.bin(tenter_coding$LocusName,tenter_coding$FST,breaks = seq(tenter_coding[1,"LocusName"],tenter_coding[nrow(tenter_coding),"LocusName"],region_size))
  # fst_tenter_coding_binned <- cbind(fst_tenter_coding_bin$stats[2,],fst_tenter_coding_bin$centers)
  # fst_tenter_coding_binned <- as.data.frame(fst_tenter_coding_binned)
  fst_coding_binned <- as.data.frame(cbind(coding$y,coding$x))
  colnames(fst_coding_binned) <- c("FST_coding","Position")
  
  # fst_tenter_noncoding_bin <- stats.bin(tenter_noncoding$LocusName,tenter_noncoding$FST,breaks = seq(tenter_noncoding[1,"LocusName"],tenter_noncoding[nrow(tenter_noncoding),"LocusName"],region_size))
  # fst_tenter_noncoding_binned <- cbind(fst_tenter_noncoding_bin$stats[2,],fst_tenter_noncoding_bin$centers)
  # fst_tenter_noncoding_binned <- as.data.frame(fst_tenter_noncoding_binned)
  fst_noncoding_binned <- as.data.frame(cbind(noncoding$y,noncoding$x))
  colnames(fst_noncoding_binned) <- c("FST_noncoding","Position")
  
  # fst_severn_coding_bin <- stats.bin(severn_coding$LocusName,severn_coding$FST,breaks = seq(severn_coding[1,"LocusName"],severn_coding[nrow(severn_coding),"LocusName"],region_size))
  # fst_severn_coding_binned <- cbind(fst_severn_coding_bin$stats[2,],fst_severn_coding_bin$centers)
  # fst_severn_coding_binned <- as.data.frame(fst_severn_coding_binned)
  # fst_severn_coding_binned <- as.data.frame(cbind(severn_coding$y,severn_coding$x))
  # colnames(fst_severn_coding_binned) <- c("FST_severn_coding","Position")
  
  # fst_severn_noncoding_bin <- stats.bin(severn_noncoding$LocusName,severn_noncoding$FST,breaks = seq(severn_noncoding[1,"LocusName"],severn_noncoding[nrow(severn_noncoding),"LocusName"],region_size))
  # fst_severn_noncoding_binned <- cbind(fst_severn_noncoding_bin$stats[2,],fst_severn_noncoding_bin$centers)
  # fst_severn_noncoding_binned <- as.data.frame(fst_severn_noncoding_binned)
  # fst_severn_noncoding_binned <- as.data.frame(cbind(severn_noncoding$y,severn_noncoding$x))
  # colnames(fst_severn_noncoding_binned) <- c("FST_severn_noncoding","Position")
  
  fst_merge <- cbind(fst_coding_binned$Position, 
                     fst_coding_binned$FST_coding,
                     fst_noncoding_binned$FST_noncoding
  )
  
  fst_merge <- as.data.frame(fst_merge)
  colnames(fst_merge) <- c("position",
                           "FST_coding",
                           "FST_noncoding"
  )
  fst_merge <- fst_merge[complete.cases(fst_merge),]
  
  #   a) Adaptive loci diverged, neutral loci did not diverged; 
  # - FST coding Severn (mean region 100Kbp) > 1 SD + FST coding Tenterfield  (mean region 100Kbp)
  # - FST coding Severn (mean region 100Kbp) > 1 SD + FST non-coding Severn  (mean region 100Kbp)
  
  scenario_a <- fst_merge[which(fst_merge$FST_coding > 
                                  (n_sd*coding_sd + fst_merge$FST_coding) &
                                  fst_merge$FST_coding > 
                                  (n_sd*noncoding_sd + fst_merge$FST_noncoding)),]
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
  
  scenario_b <- fst_merge[which(fst_merge$FST_noncoding > 
                                  (n_sd*noncoding_sd + fst_merge$FST_noncoding) &
                                  fst_merge$FST_noncoding > 
                                  (n_sd*coding_sd + fst_merge$FST_coding)),]
  
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
  
  scenario_c <- fst_merge[which(fst_merge$FST_coding > 
                                  (n_sd*coding_sd + fst_merge$FST_coding) &
                                  fst_merge$FST_noncoding > 
                                  (n_sd*noncoding_sd + fst_merge$FST_noncoding)),]
  
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
  
  for(x in 1:nrow(scenario_tmp)){
    genes_tmp <- which(scenario_tmp[x,"position"] %between% y_mRNA)
    if(length(genes_tmp) >0){
      gene_names_tmp <- mRNA[genes_tmp,]
      gene_names_tmp$chr <- chrom_n
      gene_names <- bind_rows(gene_names,gene_names_tmp)
    }else{
      next()
    }
    
  }
  
  if(nrow(scenario_tmp)>0){

    join_df_melt$stat <- factor(join_df_melt$variable, levels=c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3"))

    # cols <- c("Coding Regions"= "dodgerblue","Non-coding Regions"="dodgerblue4","Coding Regions_Tenter"="gold","Non-coding Regions_Tenter"="gold4")

    cols <- c("Coding Regions"= "lightgreen","Non-coding Regions"="darkgreen")
    #
    # cols <- c("Coding Regions"= "deeppink","Non-coding Regions"="deepskyblue3")

    f_labels <- data.frame(stat = factor("FST",levels = c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3")),
                           x = scenario_tmp$position ,
                           label =  scenario_tmp$scenario)

        p1 <-
          # ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene,fill=gene))+
          ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene))+
          facet_wrap(~stat,scales="free_y",ncol=1)+
          geom_segment(data=regions,aes(y=y, yend = y, x = start, xend = end, colour=region),linewidth = 5,  inherit.aes = FALSE) +
          stat_smooth(span = span/ span_2,se=FALSE,method = "loess",n=round(tail(join_df_melt,1)$LocusName/region_size,0))+
          theme_linedraw() +
          scale_x_continuous(breaks=seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000),labels=as.character(round(seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000)/1000000))) +
          xlab("Genome position (Mbp)")+
          labs(title = paste("Chromosome",chrom,nrow(join_df),"SNPs"), color = "") +
        scale_color_manual(values= cols) +
          theme(legend.position = "bottom",
                legend.title=element_blank(),
                axis.title.y=element_blank(),
                plot.title = element_text(hjust = 0.5)) +
          geom_label(data=f_labels,aes(label = label,fill=label),x =f_labels$x, y = -0.05,inherit.aes = FALSE,size =3,show.legend=FALSE)

         ggsave(paste0("v4_platypus_ox_chr_",chrom,".pdf"),  width = 14, height =14, units = "in", dpi="retina", bg = "transparent" )

  }else{

    join_df_melt$stat <- factor(join_df_melt$variable, levels=c("pvaluesRightTail","FST","He","PCA_Axis1","PCA_Axis2","PCA_Axis3"))

    # cols <- c("Coding Regions"= "dodgerblue","Non-coding Regions"="dodgerblue4","Coding Regions_Tenter"="gold","Non-coding Regions_Tenter"="gold4")

    cols <- c("Coding Regions"= "lightgreen","Non-coding Regions"="darkgreen")
    #
    # cols <- c("Coding Regions"= "deeppink","Non-coding Regions"="deepskyblue3")


    p1 <-
      # ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene,fill=gene))+
      ggplot(join_df_melt,aes(x=LocusName,y=value,color=gene))+
      facet_wrap(~stat,scales="free_y",ncol=1)+
      stat_smooth(span = span/ span_2,se=FALSE,method = "loess",n=round(tail(join_df_melt,1)$LocusName/region_size,0))+
      theme_linedraw() +
      scale_x_continuous(breaks=seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000),labels=as.character(round(seq(t1_order[1,"LocusName"], t1_order[nrow(t1_order),"LocusName"], 5000000)/1000000))) +
      xlab("Genome position (Mbp)")+
      labs(title = paste("Chromosome",chrom,nrow(join_df),"SNPs"), color = "") +
      scale_color_manual(values= cols) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5))

    ggsave(paste0("v4_platypus_ox_chr_",chrom,".pdf"),  width = 14, height =14, units = "in", dpi="retina", bg = "transparent" )

  }

}

write.csv(gene_names,file = "list_genes_ox.csv")
write.csv(scenario_fin,file = "scenario_fin_ox.csv")
print(n_snps)

# scenario_a = 104 - 173 - 115
# scenario_b = 44 - 82 - 46 
# scenario_c = 193 - 47 - 185 



# write.csv(list_genes,file = "list_genes.csv")
