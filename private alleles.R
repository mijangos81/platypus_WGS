library(dartR)
library(fields)
library(ggplot2)
library(patchwork)
library(ape)
library(stringr)
library(dplyr)
library(tictoc)
library(data.table)

adjust_density <- 1/2
region_size <- 200000
span <- 5000

chrom <- "1a"
gl_26 <- gl.load("gl_26_chr_1a.rds")
gl_ox <- gl.load("gl_ox_chr_1a.rds")
pa_pairwise_26 <- readRDS(file = "pa_pairwise_26_chr1a.rds")
pa_one2rest_26 <- readRDS(file = "pa_one2rest_26_chr1a.rds")
pa_pairwise_ox <- readRDS(file = "pa_pairwise_ox_chr1a.rds")
pa_one2rest_ox <- readRDS(file = "pa_one2rest_ox_chr1a.rds")

gl_ox_2 <-
  gl.filter.locmetric(
    gl_ox,
    metric = "RFGQ_ALL",
    lower = 20,
    upper = max(gl_ox$other$loc.metrics$RFGQ_ALL),
    verbose = 0
  )

tas_nq <- gl.keep.pop(gl_ox_2,pop.list = c("NQLD","TAS") )
tas_nq <- gl.filter.monomorphs(tas_nq)
fst_res <- utils.basic.stats(tas_nq)

fd_tas_nq <- as.numeric(pa_one2rest_ox$names_loci$TAS_Res$fd)
pa_nq_tas <- as.numeric(pa_pairwise_ox$names_loci$NQLD_TAS$pop1_pop2_pa)
pa_tas_nq <- as.numeric(pa_pairwise_ox$names_loci$NQLD_TAS$pop2_pop1_pa)

chrom_n <- "1a"
chr <- "1"
chrom_dir <- paste0("/Users/s441489/final_platy/chr_1/")

########################GFF########################

gff <- read.gff(paste0(chrom_dir,"chr_", chr, "_mOrnAna1.pri.v4.gff"))

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

last_snp_cod <- tail(gl_ox$position,1)
last_snp_cod <- which(non_coding$end < last_snp_cod)
last_snp_cod <- tail(last_snp_cod,1)
non_coding <- non_coding[1:last_snp_cod,]
non_coding$length <- non_coding$end - non_coding$start
#check whether is true that not real non-coding regions have negative value
non_coding <- non_coding[which(non_coding$length>0),]

last_snp_non <- tail(gl_ox$position,1)
last_snp_non <- which(coding$end < last_snp_non)
last_snp_non <- tail(last_snp_non,1)
coding <- coding[1:last_snp_non,]
coding$y <- -Inf
coding$length <- coding$end - coding$start

regions <- rbind(coding,non_coding)

fd_genes <- list(genes$attributes)
for(g in 1:nrow(genes)){
  fd_genes[[g]]<-which(fd_tas_nq %inrange% c(genes[g,"start"],genes[g,"end"]))
  
}

fd_genes_df <- data.frame( fd=unlist(lapply(fd_genes,length)),name= genes$attributes,length=genes$length,pos=genes$start)

fd_genes_df$expected <- fd_genes_df$length * exp_fd_noncoding

fd_genes_df$ratio <- fd_genes_df$fd/ fd_genes_df$expected 

fd_genes_df <- fd_genes_df[order(fd_genes_df$ratio,decreasing = T),]

###########################################################################

p1p2 <- data.frame(fd=pa_nq_tas,gene=NA)
p1p2$gene <- p1p2$fd %inrange% y_genes
p1p2[p1p2$gene==FALSE,"gene"] <- "Non-coding Regions" 
p1p2[p1p2$gene==TRUE,"gene"] <- "Coding Regions"

cols <- c("Coding Regions"= "orange","Non-coding Regions"="orange4")

p1 <- ggplot(p1p2,aes(x=fd,fill=gene,color=gene)) +
  geom_density(adjust=adjust_density, alpha=1/3) +
  # facet_wrap(~gene,scales="free_y",ncol=1)+
  theme_linedraw() +
  scale_x_continuous(breaks=seq(p1p2[1,"fd"], p1p2[nrow(p1p2),"fd"], 5000000),labels=as.character(round(seq(p1p2[1,"fd"], p1p2[nrow(p1p2),"fd"], 5000000)/1000000))) +
  xlab("Genome position (Mbp)")+
  ylab("Density private alleles")+
   labs(title = "Private alelles between North Queensland and Tasmania", color = "") +
   scale_color_manual(values= cols) +
  scale_fill_manual(values= cols) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
###########################################################################


p2p1 <- data.frame(fd=pa_tas_nq,gene=NA)
p2p1$gene <- p2p1$fd %inrange% y_genes
p2p1[p2p1$gene==FALSE,"gene"] <- "Non-coding Regions" 
p2p1[p2p1$gene==TRUE,"gene"] <- "Coding Regions"

cols <- c("Coding Regions"= "orange","Non-coding Regions"="orange4")

p2 <- ggplot(p2p1,aes(x=fd,fill=gene,color=gene)) +
  geom_density(adjust=adjust_density, alpha=1/3) +
  # facet_wrap(~gene,scales="free_y",ncol=1)+
  theme_linedraw() +
  scale_x_continuous(breaks=seq(p2p1[1,"fd"], p2p1[nrow(p2p1),"fd"], 5000000),labels=as.character(round(seq(p2p1[1,"fd"], p2p1[nrow(p2p1),"fd"], 5000000)/1000000))) +
  xlab("Genome position (Mbp)")+
  ylab("Density private alleles")+
  labs(title = "Private alelles between Tasmania and North Queensland", color = "") +
  scale_color_manual(values= cols) +
  scale_fill_manual(values= cols) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
###########################################################################


fd_df <- data.frame(fd=fd_tas_nq,gene=NA)
fd_df$gene <- fd_df$fd %inrange% y_genes
fd_df[fd_df$gene==FALSE,"gene"] <- "Non-coding Regions" 
fd_df[fd_df$gene==TRUE,"gene"] <- "Coding Regions"

fd_coding <- fd_df[which(fd_df$gene=="Coding Regions"),]
fd_noncoding <- fd_df[which(fd_df$gene=="Non-coding Regions"),]

exp_fd_coding <- nrow(fd_coding)/(sum(genes$length))
exp_fd_noncoding <- nrow(fd_noncoding)/(sum(non_coding$length))

cols <- c("Coding Regions"= "orange","Non-coding Regions"="orange4")

p3 <- ggplot(fd_df,aes(x=fd,fill=gene,color=gene)) +
  geom_density(adjust=adjust_density, alpha=1/3) +
   geom_segment(data=regions,aes(y=y, yend = y, x = start, xend = end, colour=region),linewidth = 5,  inherit.aes = FALSE) +
  # facet_wrap(~gene,scales="free_y",ncol=1)+
  theme_linedraw() +
  scale_x_continuous(breaks=seq(fd_df[1,"fd"], fd_df[nrow(fd_df),"fd"], 5000000),labels=as.character(round(seq(fd_df[1,"fd"], fd_df[nrow(fd_df),"fd"], 5000000)/1000000))) +
  xlab("Genome position (Mbp)")+
  ylab("Density fixed differences")+
  labs(title = "Fixed differences between Tasmania and North Queensland", color = "") +
  scale_color_manual(values= cols) +
   scale_fill_manual(values= cols) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

###########################################################################


stats_tas_nq <- data.frame(pos=tas_nq$position,fst=fst_res$perloc$Fst,he=fst_res$perloc$Hs)

stats_tas_nq <- stats_tas_nq[complete.cases(stats_tas_nq$fst),]
stats_tas_nq$gene <- stats_tas_nq$pos %inrange% y_genes
stats_tas_nq[stats_tas_nq$gene==FALSE,"gene"] <- "Non-coding Regions" 
stats_tas_nq[stats_tas_nq$gene==TRUE,"gene"] <- "Coding Regions"

join_df_melt <- melt(stats_tas_nq,id.vars=c("pos","gene"))

span_2 <-  nrow(stats_tas_nq)
  

  p4 <-
    ggplot(join_df_melt,aes(x=pos,y=value,color=gene))+
    facet_wrap(~variable,scales="free_y",ncol=1)+
    stat_smooth(span = span/span_2,se=FALSE,method = "loess",n=round(tail(join_df_melt,1)$pos/region_size,0))+
    theme_linedraw() +
     scale_x_continuous(breaks=seq(stats_tas_nq[1,"pos"], stats_tas_nq[nrow(stats_tas_nq),"pos"], 5000000),labels=as.character(round(seq(stats_tas_nq[1,"pos"], stats_tas_nq[nrow(stats_tas_nq),"pos"], 5000000)/1000000))) +
    xlab("Genome position (Mbp)")+
    ylab("Value")+
    # labs(title = paste("Chromosome",chrom,nrow(join_df),"SNPs"), color = "") +
    scale_color_manual(values= cols) +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  print(
    p3/p4 + plot_layout(ncol =1 , heights = c(1,2))
  )

  ggsave(paste0("pa_ox_chr_",chrom,".pdf"),  width = 12, height =12, units = "in", dpi="retina", bg = "transparent" )








