
# p <- colMeans(as.matrix(x), na.rm = TRUE) / 2
# # ignore loci with just missing data
# i0 <- which(!is.na(p))  
# logp <- ifelse(!is.finite(log(p)), 0, log(p))
# log1_p <- ifelse(!is.finite(log(1 - p)), 0, log(1 - p))
# one_H_alpha_all <- -(p * logp + (1 - p) * log1_p)
# 
# x2 <- seppop(x)
# 
# 
#   i1 <- which(!is.na(colMeans(as.matrix(x2[[1]]), na.rm = TRUE) / 2))
#   i2 <- which(!is.na(colMeans(as.matrix(x2[[2]]), na.rm = TRUE) / 2))
#   tt <- table(c(i0, i1, i2))
#   index <- as.numeric(names(tt)[tt == 3])
#   dummys <-
#     one_H_alpha_all[i0 %in% index] - 
#     (one_H_alpha_es[[x[1]]]$dummys[i1 %in% index] + 
#        one_H_alpha_es[[x[2]]]$dummys[i2 %in% index]) / 2


mutual_information <- function(mat) {
  EstMLEFun <- function(mat) {
    # MLE
    entropyFun <- function(p) {
      p <- p[p > 0]
      out <- -sum(p * log(p))
      return(out)
    }
    n <- sum(mat)
    prob.hat <- mat / n
    px.hat <- apply(prob.hat, 1, sum)
    py.hat <- apply(prob.hat, 2, sum)
    I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
    # MLE of Mutual Information!
    return(I.hat)
  }
  mydata <- as.matrix(mat)
  est <- EstMLEFun(mydata)
  return(est)
}

x <- platypus.gl
x <- gl.drop.pop(x,pop.list = "TENTERFIELD")
x <- gl.filter.callrate(x,threshold = 1)
x <- gl.filter.monomorphs(x)

tmp_genind <- dartR::gl2gi(tmp)
tmp_hierf <-  hierfstat::genind2hierfstat(tmp_genind)
tmp_allele_matrix <- allele.count(tmp_hierf)
tmp_MI <- lapply(tmp_allele_matrix, mutual_information)
res[i, "MI"] <- round(mean(unlist(tmp_MI)), 3)
res[i, "MI_SE"] <- round(std(unlist(tmp_MI)), 3)

res <- mutual_information(x)
res_2 <- gl.report.diversity(x)

library(hierfstat)

tmp_genind <- dartR::gl2gi(x)
tmp_hierf <-  hierfstat::genind2hierfstat(tmp_genind)
tmp_allele_matrix <- allele.count(tmp_hierf)
tmp_MI <- lapply(tmp_allele_matrix, mutual_information)
res <- mean(round(unname(unlist(tmp_MI)),3 ))

res_2 <- gl.report.diversity(x)

# [1] 0.01885066


shannon <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  out <- - sum(p * log(p))
  return(out)
}

MI <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  
  mat <- as.matrix(x)
  
  n <- sum(mat)
  prob.hat <- mat / n
  px.hat <- apply(prob.hat, 1, sum)
  py.hat <- apply(prob.hat, 2, sum)
  I.hat <- shannon(px.hat) + shannon(py.hat) - shannon(prob.hat)
  # MLE of Mutual Information!
  return(I.hat)
}

AFD <- function(pop1,pop2){
  
  p1 <- as.matrix(pop1)
  p2 <- as.matrix(pop2)
  p1alf <- colMeans(p1, na.rm = TRUE) / 2
  p2alf <- colMeans(p2, na.rm = TRUE) / 2
  
  return(abs(p1alf - p2alf))
}

AA <- function(pop1,pop2){
  
  p1 <- as.matrix(pop1)
  p2 <- as.matrix(pop2)
  p1alf <- colMeans(p1, na.rm = TRUE) / 2
  p2alf <- colMeans(p2, na.rm = TRUE) / 2
  
  pmax <- apply(cbind(p1alf,p2alf),1,max)
  BCAFD <- abs(p1alf - p2alf)
  AA_tmp <- BCAFD / (0.6152 + (0.3985 * pmax))
  
  return(AA_tmp)
}


FST <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Fst <- utils.basic.stats(x)
  
  return(Fst$perloc$Fst)
  
}

FSTp <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Fstp <- utils.basic.stats(x)
  
  return(Fstp$perloc$Fstp)
  
}

Dest <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Dest_tmp <- utils.basic.stats(x)
  
  return(Dest_tmp$perloc$Dest)
  
}

Gst_H <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Gst_H_tmp <- utils.basic.stats(x)
  
  return(Gst_H_tmp$perloc$Gst_H)
  
}








