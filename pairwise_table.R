

#' @param x A genlight object containing the SNP genotypes [required].
#' @param pw_fun A function whose inputs are pop1 and pop2[required].
#' @param dec Number of decimals [default 4].

pairwise_table <- function(x,
                           pw_fun,
                           dec = 4){
  
  pop_names <- popNames(x)
  x <- seppop(x)
  ix <- setNames(seq_along(pop_names), pop_names)
  pp <- outer(ix[-1L], ix[-length(ix)], 
              function(ivec, jvec){
                vapply(seq_along(ivec), 
                       function(k) {
                         i <- ivec[k]
                         j <- jvec[k]
                         if (i > j){
                           # x_tmp <- rbind(x[[i]],x[[j]])
                           round(mean(pw_fun(x[[i]],x[[j]]), na.rm = TRUE), dec)
                         }else{
                           NA_real_
                         }
                       }, numeric(1))
              })
  return(pp)
  
}