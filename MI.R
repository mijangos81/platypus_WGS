
#Convert data set into matrix and calculate column means divided by 2. 
#Remove missing values and store result as new vector 'p'.
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

#Explanation of the function:

#The outer function ‘mutual_information’ takes a single argument mat, which is the joint probability distribution 
#matrix.
#Inside the ‘mutual_information’ function, there is a nested function called ‘EstMLEFun’.

#This function is used to calculate the mutual information using the MLE approach.

#‘entropyFun’ is a helper function defined inside ‘EstMLEFun’, used to calculate the entropy of a probability 
#distribution ‘p’. Entropy is a measure of uncertainty or randomness in a probability distribution. 

#It's defined as ‘-sum(p * log(p))’, where the sum is taken over non-zero elements of the probability distribution.

#‘n’ is the total count (sum) of all elements in the joint probability distribution matrix mat.
#‘prob.hat’ calculates the estimated probability distribution matrix by dividing each element of ‘mat’ by ‘n’.

#‘px.hat’ and ‘py.hat’ calculate the marginal probability distributions for the rows and columns of ‘prob.hat’,
#respectively, using the ‘apply()’ function. The ‘apply()’ function is used to apply a function 
# (in this case, the ‘sum()’ function) over rows or columns of a matrix.

#‘I.hat’ calculates the mutual information estimate using MLE. It adds the entropies of the row distribution 
#‘px.hat’ and the column distribution ‘py.hat’, and then subtracts the entropy of the joint distribution ‘prob.hat’.

#The MLE estimate of mutual information, ‘I.hat’, is returned by the ‘EstMLEFun’ function.

#The main function ‘mutual_information’ converts the input mat to a matrix (if it's not already) using ‘as.matrix()’, 
#then calls the nested function ‘EstMLEFun’ with the matrix data.

#The mutual information estimate (‘est’) obtained from ‘EstMLEFun’ is returned as the final result of the ‘mutual_information’ function.

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

# x <- platypus.gl
# x <- gl.drop.pop(x,pop.list = "TENTERFIELD")
# x <- gl.filter.callrate(x,threshold = 1)
# x <- gl.filter.monomorphs(x)
# 
# tmp_genind <- dartR::gl2gi(x)
# tmp_hierf <-  hierfstat::genind2hierfstat(tmp_genind)
# tmp_allele_matrix <- allele.count(tmp_hierf)
# tmp_MI <- lapply(tmp_allele_matrix, mutual_information)
# res[i, "MI"] <- round(mean(unlist(tmp_MI)), 3)
# res[i, "MI_SE"] <- round(std(unlist(tmp_MI)), 3)
# 
# res <- mutual_information(x)
# res_2 <- gl.report.diversity(x)

library(hierfstat)
# 
# tmp_genind <- dartR::gl2gi(x)
# tmp_hierf <-  hierfstat::genind2hierfstat(tmp_genind)
# tmp_allele_matrix <- allele.count(tmp_hierf)
# tmp_MI <- lapply(tmp_allele_matrix, mutual_information)
# res <- mean(round(unname(unlist(tmp_MI)),3 ))
# 
# res_2 <- gl.report.diversity(x)

# [1] 0.01885066


#Explanation of the code:

#‘shannon <- function(x) {‘: This line defines a function named ‘shannon’ that takes one argument ‘x’. 
#The argument ‘x’ is expected to be a numeric vector representing a probability distribution (e.g., probabilities of events).

#‘x <- x[x > 0]’: This line filters the input vector ‘x’ to remove any elements that are less than or equal to zero. 
#The purpose of this step is to exclude any zero probabilities, as the Shannon entropy formula includes a logarithm operation, 
#which is undefined for zero values.

#‘p <- x / sum(x)’: This line calculates the normalized probabilities (proportions) of the elements in the filtered vector ‘x’. 
#Each element in ‘x’ is divided by the sum of all elements in ‘x’ to ensure that the probabilities sum up to 1.

#‘out <- -sum(p * log(p))’: This line calculates the Shannon entropy. It multiplies each probability ‘p’ with the logarithm 
#(base 2 by default) of that probability and then takes the negative sum of these values. 
#The negative sign is used to ensure that the entropy value is positive or zero.

#‘return(out)’: This line returns the calculated Shannon entropy value out as the result of the function.
    
shannon <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  out <- - sum(p * log(p))
  return(out)
}

#Explanation of the code:

#‘MI <- function(pop1, pop2) {‘: This line defines a function named ‘MI’ that takes two arguments, ‘pop1’ and ‘pop2’, 
#representing the two populations for which mutual information needs to be estimated.

#‘x <- rbind(pop1, pop2)’: This line creates a new matrix ‘x’ by binding the two populations ‘pop1’ and ‘pop2’ 
#together as rows. The resulting matrix ‘x’ will have two rows, where each row represents a population.

#‘mat <- as.matrix(x)’: This line converts the matrix ‘x’ into a regular matrix format using the ‘as.matrix()’ function. 
#Converting ‘x’ into a matrix is important for the subsequent calculations.

#‘n <- sum(mat)’: This line calculates the total count (sum) of all elements in the joint probability distribution 
#matrix ‘mat’

#‘prob.hat <- mat / n’: This line calculates the estimated joint probability distribution matrix by dividing each 
#element of ‘mat’ by ‘n’. The resulting matrix ‘prob.hat’ contains the estimated probabilities of joint events.

#‘px.hat <- apply(prob.hat, 1, sum)’: This line calculates the marginal probability distribution for ‘pop1’ by 
#summing the probabilities along rows of ‘prob.hat’ using the ‘apply()’ function. The result is a vector containing 
#the estimated probabilities for ‘pop1’.

#‘py.hat <- apply(prob.hat, 2, sum)’: This line calculates the marginal probability distribution for ‘pop2’ by 
#summing the probabilities along columns of ‘prob.hat’ using the ‘apply()’ function. The result is a vector containing 
#the estimated probabilities for ‘pop2’.

#‘I.hat <- shannon(px.hat) + shannon(py.hat) - shannon(prob.hat)’: This line calculates the mutual information 
#estimate using the Shannon entropy (‘shannon’) function. It adds the entropies of the marginal distributions ‘px.hat’ 
#and ‘py.hat’, and then subtracts the entropy of the joint distribution ‘prob.hat’.

#The MLE estimate of mutual information,’ I.hat’, is returned as the result of the MI function.
    

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

#Explanation of the code:
  
#‘AFD <- function(pop1, pop2) {‘: This line defines a function named ‘AFD’ that takes two arguments, ‘pop1 and pop2’, 
#representing the two populations for which the Absolute Frequency Difference needs to be calculated.

#‘p1 <- as.matrix(pop1)’: This line converts the input ‘pop1’ into a matrix format using the ‘as.matrix()’ function. 
#Converting ‘pop1’ into a matrix is necessary for subsequent calculations.

#‘p2 <- as.matrix(pop2)’: This line converts the input ‘pop2’ into a matrix format using the ‘as.matrix()’ function. 
#Converting ‘pop2’ into a matrix is necessary for subsequent calculations.

#‘p1alf <- colMeans(p1, na.rm = TRUE) / 2’: This line calculates the average frequencies for each column in ‘pop1’ and 
#then divides them by 2. The ‘colMeans()’ function is used to compute the mean value for each column of the matrix ‘p1’,
#and the ‘na.rm = TRUE’ argument ensures that any missing values (NA) are ignored in the calculation.

#‘p2alf <- colMeans(p2, na.rm = TRUE) / 2’: This line calculates the average frequencies for each column in ‘pop2’ and 
#then divides them by 2. The ‘colMeans()’ function is used to compute the mean value for each column of the matrix ‘p2’,
#and the ‘na.rm = TRUE’ argument ensures that any missing values (NA) are ignored in the calculation.

#‘return(abs(p1alf - p2alf))’: This line returns the absolute difference between the average frequencies calculated 
#for ‘pop1’ and ‘pop2’. The ‘abs()’ function ensures that the result is a positive value.
    

AFD <- function(pop1,pop2){
  
  p1 <- as.matrix(pop1)
  p2 <- as.matrix(pop2)
  p1alf <- colMeans(p1, na.rm = TRUE) / 2
  p2alf <- colMeans(p2, na.rm = TRUE) / 2
  
  return(abs(p1alf - p2alf))
}

#Explanation of the code:

#‘AA <- function(pop1, pop2) {‘: This line defines a function named ‘AA’ that takes two arguments, ‘pop1’ and ‘pop2’, 
#representing the two populations for which the Adjusted Absolute Frequency Difference needs to be calculated.

#‘p1 <- as.matrix(pop1’): This line converts the input ‘pop1’ into a matrix format using the ‘as.matrix()’ function. 
#Converting ‘pop1’ into a matrix is necessary for subsequent calculations.

#‘p2 <- as.matrix(pop2)’: This line converts the input ‘pop2’ into a matrix format using the ‘as.matrix()’ function. 
#Converting ‘pop2’ into a matrix is necessary for subsequent calculations.

#‘p1alf <- colMeans(p1, na.rm = TRUE) / 2’: This line calculates the average frequencies for each column in ‘pop1’ and 
#then divides them by 2. The ‘colMeans()’ function is used to compute the mean value for each column of the matrix ‘p1’, and the ‘na.rm = TRUE’ argument ensures that any missing values (NA) are ignored in the calculation.

#‘p2alf <- colMeans(p2, na.rm = TRUE) / 2’: This line calculates the average frequencies for each column in ‘pop2’ and 
#then divides them by 2. The ‘colMeans()’ function is used to compute the mean value for each column of the matrix ‘p2’, and the ‘na.rm = TRUE’ argument ensures that any missing values (NA) are ignored in the calculation.

#‘pmax <- apply(cbind(p1alf, p2alf), 1, max)’: This line calculates the maximum proportion of alleles in each column 
#by taking the maximum value between the corresponding elements of ‘p1alf’ and ‘p2alf’. The ‘cbind()’ function is used to combine the vectors ‘p1alf’ and ‘p2alf’ column-wise, and the ‘apply()’ function with ‘1’ as the second argument ensures that the ‘max()’ function is applied along rows.

#‘BCAFD <- abs(p1alf - p2alf)’: This line calculates the Absolute Frequency Difference (AFD) between the populations 
#by taking the absolute difference between the average frequencies of ‘pop1’ and ‘pop2’.

#‘AA_tmp <- BCAFD / (0.6152 + (0.3985 * pmax))’: This line calculates the Adjusted Absolute Frequency Difference (AA) 
#by dividing the AFD ‘BCAFD’ by a formula that includes ‘pmax’ and constant values 0.6152 and 0.3985. The formula adjusts the AFD based on the proportion of alleles and the weighting factors.

#The calculated Adjusted Absolute Frequency Difference, ‘AA_tmp’, is returned as the result of the ‘AA’ function.
    
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

#Explanation of the code:

#‘FST <- function(pop1, pop2) {‘: This line defines a function named ‘FST’ that takes two arguments, ‘pop1’ and ‘pop2’,
#representing the two populations for which the FST needs to be calculated.
    
#‘x <- rbind(pop1, pop2)’: This line creates a new matrix ‘x’ by binding the two populations ‘pop1’ and ‘pop2’ together
#as rows. The resulting matrix ‘x’ will have two rows, where each row represents a population.
    
#‘Fst <- utils.basic.stats(x)’: This line calls the custom function ‘utils.basic.stats()’ with the matrix ‘x’ as an 
#argument. The function is assumed to calculate various basic genetic statistics, including the FST.
    
#‘return(Fst$perloc$Fst)’: This line returns the FST value calculated by the ‘utils.basic.stats()’ function. 
#The FST value is extracted from the result of the function using ‘$perloc$Fst’. The exact structure of the result 
#returned by ‘utils.basic.stats()’ is not provided, so we can only infer that the FST value is stored in ‘$perloc$Fst’ 
#based on the code.

#To use this function, you would need to have the ‘utils.basic.stats()’ function available, which is not part of the 
#standard R library. If you don't have this function or its definition, you will need to either implement it or obtain it from another source to use the FST function.


FST <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Fst <- utils.basic.stats(x)
  
  return(Fst$perloc$Fst)
  
}

#‘FSTp <- function(pop1, pop2) {‘: This line defines a function named ‘FSTp’ that takes two arguments, ‘pop1’ and ‘pop2’, 
#representing the two populations for which the FST and its p-value need to be calculated.

#‘x <- rbind(pop1, pop2)’: This line creates a new matrix ‘x’ by binding the two populations ‘pop1’ and ‘pop2’ together 
#as rows. The resulting matrix ‘x’ will have two rows, where each row represents a population.

#‘Fstp <- utils.basic.stats(x)’: This line calls the custom function ‘utils.basic.stats()’ with the matrix ‘x’ as an 
#argument. The function is assumed to calculate various basic genetic statistics, including the FST and its associated p-value.
  
#‘return(Fstp$perloc$Fstp)’: This line returns the p-value associated with the FST value calculated by the ‘utils.basic.stats()’ 
#function. The p-value is extracted from the result of the function using ‘$perloc$Fstp’. The exact structure of the result returned by ‘utils.basic.stats()’ is not provided, so we can only infer that the p-value is stored in ‘$perloc$Fstp’ based on the code.
  

FSTp <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Fstp <- utils.basic.stats(x)
  
  return(Fstp$perloc$Fstp)
  
}

#Explanation of the code:

#‘Dest <- function(pop1, pop2) {‘: This line defines a function named ‘Dest’ that takes two arguments, ‘pop1’ and ‘pop2’,
#representing the two populations for which the Dest needs to be calculated.
    
#‘x <- rbind(pop1, pop2)’: This line creates a new matrix ‘x’ by binding the two populations ‘pop1’ and ‘pop2’ together 
#as rows. The resulting matrix ‘x’ will have two rows, where each row represents a population.
    
#‘Dest_tmp <- utils.basic.stats(x)’: This line calls the custom function ‘utils.basic.stats()’ with the matrix ‘x’ as 
#an argument. The function is assumed to calculate various basic genetic statistics, including the Dest.
    
#‘return(Dest_tmp$perloc$Dest)’: This line returns the Dest value calculated by the ‘utils.basic.stats()’ function. 
#The Dest value is extracted from the result of the function using ‘$perloc$Dest’. The exact structure of the result returned by ‘utils.basic.stats()’ is not provided, so we can only infer that the Dest value is stored in ‘$perloc$Dest’ based on the code.
    
Dest <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Dest_tmp <- utils.basic.stats(x)
  
  return(Dest_tmp$perloc$Dest)
  
}

#Explanation of the code:

#‘Gst_H <- function(pop1, pop2) {‘: This line defines a function named ‘Gst_H’ that takes two arguments, ‘pop1’ and ‘pop2’, 
#representing the two populations for which the Gst_H needs to be calculated.
    
#‘x <- rbind(pop1, pop2)’: This line creates a new matrix ‘x’ by binding the two populations ‘pop1’ and ‘pop2’ together 
#as rows. The resulting matrix ‘x’ will have two rows, where each row represents a population.
    
#‘Gst_H_tmp <- utils.basic.stats(x)’: This line calls the custom function ‘utils.basic.stats()’ with the matrix ‘x’ as 
#an argument. The function is assumed to calculate various basic genetic statistics, including Gst_H.
    
#‘return(Gst_H_tmp$perloc$Gst_H)’: This line returns the Gst_H value calculated by the ‘utils.basic.stats()’ function. 
#The Gst_H value is extracted from the result of the function using ‘$perloc$Gst_H’. The exact structure of the result returned by ‘utils.basic.stats()’ is not provided, so we can only infer that the Gst_H value is stored in ‘$perloc$Gst_H’ based on the code.
  

Gst_H <- function(pop1,pop2){
  
  x <- rbind(pop1,pop2)  
  Gst_H_tmp <- utils.basic.stats(x)
  
  return(Gst_H_tmp$perloc$Gst_H)
  
}








