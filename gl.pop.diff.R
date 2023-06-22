#' @name gl.pop.diff
#' @title 
#' @description
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param nboots [default 0].
#' @param parallel The type of parallel operation to be used (if any). If missing, the
#'  default is taken from the option "boot.parallel" (and if that is not set, "no").
#' @param ncpus	number of processes to be used in parallel operation: typically 
#' one would chose this to the number of available CPUs.
#' @param ... 
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors A color palette or a list with as many colors as there are 
#' populations in the dataset [default discrete_palette].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details 
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife. 
#' Annals of Statistics 7, 1–26.
#' Efron [1979] suggested that Bootstraps should be at least 200.)
#' @return
#' @author 
#' @examples
#' @family report functions
#' @references
#' @export

gl.pop.diff <- function(x,
                        nboots = 0,
                        conf = 0.95,
                        parallel = "no",
                        ncpus = 1,
                        digits = 4,
                        plot.out = TRUE,
                        plot_theme = theme_dartR(),
                        plot_colors = discrete_palette,
                        save2tmp = FALSE,
                        verbose = NULL
                        ){
  
  # bootstrapping function
  pop.diff <- function(x,indices){
    
    x2 <- x[,indices]
    pop_diff <- utils.basic.stats(x2,digits=digits)
    
    return(pop_diff$overall)
    
  }
  
  pops <- seppop(x)
  npops <- length(pops)
  pairs_pops <- t(combn(npops, 2))
  
  ### pairwise 
  if(npops>2){
    
    pairpop <- apply(pairs_pops, 1, function(y) {
      
      tpop <- rbind(pops[[y[1]]],pops[[y[2]]])
      
      if(nboots>0){
        
        res <- boot(data = tpop, statistic = pop.diff, R = nboots)
        
      }else{
        
        res <- utils.basic.stats(tpop,digits=digits)$overall
      }
      
      return(res)
      
    })
    
  }else{
    
      tpop <- rbind(pops[[1]],pops[[2]])
      
      if(nboots>0){
        
        pairpop <- boot(data = tpop, statistic = pop.diff, R = nboots)
        
      }else{
        
        pairpop <- utils.basic.stats(tpop,digits=digits)$overall
      }
  }

  if(nboots>0){
    
    mat_pops <- rep(list(matrix(NA, nrow = npops, ncol = npops)),12)
    
    # if more than 2 pops 
    if(npops>2){
      pairpop_2 <- Reduce(cbind,lapply(pairpop, "[[",1))
    }else{
      pairpop_2 <- as.data.frame(pairpop[[1]])
    }
    
    for(i in 1:length(mat_pops)){
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop_2[i,]
      colnames(mat_pops[[i]]) <- rownames(mat_pops[[i]]) <- names(pops)
    }
    
    names(mat_pops) <- rownames(pairpop_2)
    
    # calculation of confident intervals are based on the function Bootstrap.CI 
    # from Anne Chao's package spadeR 
    if(npops>2){
      pairpop_boot <- lapply(pairpop,"[[","t")
      pairpop.mean <- lapply(pairpop_boot,apply,2,mean)
      LCI <- as.list(1:length(pairpop.mean))
      UCI <- as.list(1:length(pairpop.mean))
    }else{
      pairpop_boot <- pairpop$t
      pairpop.mean <- colMeans(pairpop_boot)
      LCI <- NULL
      UCI <- NULL
    }

      for (i in 1:length(pairpop.mean)) {
        if(npops>2){
          LCI_tmp <- -apply(pairpop_boot[[i]],2,function(x)quantile(x,probs = (1-conf)/2)) + pairpop.mean[[i]]
          UCI_tmp <-  apply(pairpop_boot[[i]],2,function(x)quantile(x,probs = 1-(1-conf)/2)) - pairpop.mean[[i]]
          LCI[[i]] <- pairpop_2[,i] - LCI_tmp
          UCI[[i]] <- pairpop_2[,i] + UCI_tmp
        }else{
          # LCI_tmp <- -apply(pairpop_boot,2,function(x)quantile(x,probs = (1-conf)/2)) + pairpop.mean[i]
          # UCI_tmp <-  apply(pairpop_boot,2,function(x)quantile(x,probs = 1-(1-conf)/2)) - pairpop.mean[i]
          # LCI <- pairpop_2 - LCI_tmp
          # UCI <- pairpop_2 + UCI_tmp
          LCI_tmp <- apply(pairpop_boot,2,function(x)quantile(x,probs = (1-conf)/2))
          UCI_tmp <- apply(pairpop_boot,2,function(x)quantile(x,probs = 1-(1-conf)/2))
          LCI <- LCI_tmp
          UCI <- UCI_tmp
        }
      }
      
    if(npops>2){
      stat_pop <- lapply(pairpop, "[[","t0")
      CI <- Map(cbind,stat_pop,LCI,UCI)
    }else{
      stat_pop <- pairpop$t0
      CI <- Reduce(cbind,list(stat_pop,LCI,UCI))
    }
      
  }else{
    
    if(npops>2){
      mat_pops <- rep(list(matrix(NA, nrow = npops, ncol = npops)),nrow(pairpop))
    }else{
      mat_pops <- rep(list(matrix(NA, nrow = npops, ncol = npops)),length(pairpop))
    }
    
    for(i in 1:length(mat_pops)){
      mat_pops[[i]][lower.tri(mat_pops[[i]])] <- pairpop[i]
      colnames(mat_pops[[i]]) <- rownames(mat_pops[[i]]) <- names(pops)
    }

    if(npops>2){
      names(mat_pops) <- rownames(pairpop)
    }else{
      names(mat_pops) <- names(pairpop)
    }
    
  }
  
  if(nboots>0){
    return(list(mat_pops,CI))
  }else{
    return(mat_pops)
  }

  }


  
