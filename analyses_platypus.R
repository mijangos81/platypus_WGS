# loading packages
library(ape)
library(dartR)
library(stringr)
library(dplyr)
library(tictoc)
library(data.table)

########################################################################
###########################  loading datasets ##########################
########################################################################
# Objects loaded with function load_platy are: 
# genes 
# gff 
# gl_26 
# gl_ox  
# gl_dart  
# gl_dart_26  
source_load <- paste0(getwd(),"/scripts/load_platy.R")
source(source_load)
# the only parameter of the function is the chromosome name
load_list <- load_platy("X1a")
########################################################################
###########################  analyses ##########################
########################################################################
gff <- load_list$gff
genes <- load_list$genes
gene_ox <- load_list$gl_ox
gene_dart <- load_list$gl_dart
gene_26 <- load_list$gl_26
gene_dart_26 <- load_list$gl_dart_26

MHC_name <- "MHC"
mhc <- gff[gff$type=="gene",]
mhc <- mhc[grep(MHC_name,mhc$attributes),]
mhc <- cbind(mhc, as.data.frame(str_split(mhc$attributes, pattern = ";", 
                                          simplify = TRUE)))

loci <- gene_26$position
mhc_genes <- NULL
for(i in 1:nrow(mhc)){
  genes <- mhc[i,]
  y  <- list(genes$start,genes$end)
  mhc_genes_tmp <- which(loci %inrange% y)
  mhc_genes <- c(mhc_genes,mhc_genes_tmp)
}

gene_26_mhc <- gl.keep.loc(gene_26,loc.list = locNames(gene_26)[mhc_genes])

dataset_list <- list(dart = gene_dart_26, gl_26 = gene_26_mhc)
names_datset <- names(dataset_list)

for(i in 1:length(dataset_list)){
  db <- dataset_list[[i]]
  gl.report.heterozygosity(db,verbose = 0)
  ggsave(paste0("het_mhc_",names_datset[i],".pdf"),  width = 6, height = 6,
         units = "in", dpi="retina",bg = "transparent" )
}

