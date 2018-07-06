#' A meta-analysis approach with filtering for identifying gene-level
#' gene-environment interactions with genetic association data
#' 
#' This function first conducts a meta-filtering test to filter out unpromising
#' SNPs. It then runs a test of omnibus-filtering-based GxE meta-analysis
#' (ofGEM) that combines the strengths of the fixed- and random-effects
#' meta-analysis with meta-filtering. It can also analyze data from multiple
#' ethnic groups. The p-values are calculated using a sequential sampling
#' approach.
#' 
#' 
#' @param Z a matrix of test statistics for gene-environment interactions (GxE)
#' from consortium data. Each row corresponds to a SNP in a set (e.g., a gene),
#' and each column represents a study. For multi-ethnic groups, Z is a list
#' with each element being the matrix for each ethnic group.
#' @param X a matrix of filtering statistics for GxE. Each row corresponds to a
#' SNP in a set, and each column represents a study. For multi-ethnic groups, X
#' is a list with each element being the matrix for each ethnic group.
#' @param R the correlation matrix of test statistics for SNPs in a set. One
#' may use the genotype LD matrix for the set of SNPs to approximate it. This
#' matrix is used when sampling correlated testing and filtering statistics
#' under the null hypothesis and to obtain the null meta-analysis statistics.
#' For multi-ethnic groups, R is a list with each element being the correlation
#' matrix for each ethnic group.
#' @param weight the weight vector for each study, or the weight matrix for
#' each SNP and each study. If the weight is the same across SNPs, it is a
#' vector with length equaling to the number of studies. If the weight is
#' different for different SNPs, it is a matrix with each row corresponding to
#' each SNP and each column representing each study.
#' @param threshold a fixed p-value threshold for filtering test. The default
#' is 0.1.
#' @param maxSim the maximum number of samplings performed in obtaining the
#' sets of correlated testing and filtering statistics under the null. The
#' default is 1e6. This number determines the precision of the p-value
#' calculation.
#' @param tol the tolerance number to stop the sequential sampling procedure.
#' The default is 10. We count the number of sampling-based null
#' meta-statistics that is more extreme than the observed meta-statistics. We
#' sequentially increase the number of sampling with an increment of 100. The
#' sequential sampling will stop if the cumulative count reaches tol. The idea
#' is to stop pursuing a higher precision with more sampling of null if the
#' p-value appears to be not significant. If tol = 0, the number of samplings
#' equals to maxSim.
#' @return A list containing \item{pval_random_mf}{the p-value based on the
#' random-effects meta-analysis test with its corresponding meta-filtering.}
#' \item{pval_fixed_mf}{the p-value based on the fixed-effects meta-analysis
#' test with its corresponding meta-filtering.} \item{pval_ofGEM}{the p-value
#' based on using Fisher's method to aggregating the p-values of fixed- and
#' random-effects meta-analysis tests with meta-filtering} \item{nsim}{the
#' number of samplings being performed.}
#' @references Wang, Jiebiao, Qianying Liu, Brandon L. Pierce, Dezheng Huo, 
#' Olufunmilayo I. Olopade, Habibul Ahsan, and Lin S. Chen. 
#' "A meta-analysis approach with filtering for identifying gene-level gene-environment interactions."
#'  Genetic epidemiology (2018). https://doi.org/10.1002/gepi.22115
#' @examples
#' 
#' 
#' data(sim_dat)
#' 
#' pval = ofGEM(Z = sim_dat$Z, X = sim_dat$X, R = sim_dat$R, weight = rep(1/6, 6))
#' 
#' 
#' @export ofGEM

ofGEM = function(Z, X, R, weight, threshold = 0.1, maxSim = 1e6, tol = 10) {
  
  # if multi-ethnic
  
  if(is.list(Z)) {
    
    Nsnp = nrow(X[[1]])
    
    T_fixed_mf = T_random_mf = 0
    
    for(i in 1:length(Z)) {
      
      # test if weight is a matrix first
      
      weight_mat = matrix(rep(weight[[i]], Nsnp),byrow=T,nrow=Nsnp)
      
      Z_weighted = Z[[i]] * weight_mat
      
      # MF_fixed
      
      MF_fixed = rowSums(X[[i]] * weight_mat)
      
      # MF_random
      
      MF_random = rowSums(X[[i]]^2 * weight_mat^2)
      
      ## test statistics
      
      T_fixed_mf = T_fixed_mf + sum(rowSums(Z_weighted)^2 * (abs(MF_fixed) > qnorm((1 - threshold/2))))
      
      T_random_mf = T_random_mf + sum(rowSums(Z_weighted^2) * 
                          (sapply(Nsnp, function(q) davies(MF_random[q], lambda=weight_mat[q,]^2)$Qq) < threshold))
    }
    
  } else {
    
    Nsnp = nrow(X)
    
    weight_mat = matrix(rep(weight, Nsnp),byrow=T,nrow=Nsnp)
    
    Z_weighted = Z * weight_mat
    
    # MF_fixed
    
    MF_fixed = rowSums(X * weight_mat)
    
    # MF_random
    
    MF_random = rowSums(X^2 * weight_mat^2)
    
    T_fixed_mf = sum(rowSums(Z_weighted)^2 * (abs(MF_fixed) > qnorm((1 - threshold/2))))
    
    T_random_mf = sum(rowSums(Z_weighted^2) * 
                        (sapply(Nsnp, function(q) davies(MF_random[q], lambda=weight_mat[q,]^2)$Qq) < threshold))
  }
  
  
  ## Calculate p-values with the sequential sampling strategy
  
  pval = rep(1, 2)
  nsim = 0
  
  if (!all(c(T_random_mf, T_fixed_mf) == 0, na.rm = T)) {
    
    nsim = 100
    count = sim(n = nsim, R, weight, threshold, T_random_mf, T_fixed_mf)
    
    while (!all(count >= tol, na.rm = T) & nsim < maxSim) {
      
      count = count + sim(n=100, R, weight, threshold, T_random_mf, T_fixed_mf)
      nsim = nsim + 100
      
    }
    
    pval = count / nsim
    
  }
  
  # oGEMf
  
  return(list(nsim = nsim, pval_random_mf = pval[1], pval_fixed_mf = pval[2], pval_ofGEM = Fisher.test(pval)$p.value))
  
}



# simulate the statistics under the null

sim = function(n, R, weight, threshold, T_random_mf, T_fixed_mf) {
  
  T_random_mf_null = T_fixed_mf_null = NULL
  
  for (i in 1:n) {
    
    # if multi-ethnic
    
    if(is.list(R)) {
      
      # X = Z = matrix(NA, ncol(R[[1]]), nstudy)
      
      T_fixed_mf_null[i] = T_random_mf_null[i] = 0
      
      for(j in 1:length(R)) {
        
        nstudy = length(weight[[1]])
        
        weight_mat = matrix(rep(weight[[j]], ncol(R[[j]])),byrow=T,nrow=ncol(R[[j]]))
        
        ## for study-specific R: ignore it first
        # X[,j] = t(mvrnorm(1, rep(0, ncol(R[[j]])), R[[j]]))
        # Z[,j] = t(mvrnorm(1, rep(0, ncol(R[[j]])), R[[j]]))
        
        X = t(mvrnorm(nstudy, rep(0, ncol(R[[j]])), R[[j]]))
        Z = t(mvrnorm(nstudy, rep(0, ncol(R[[j]])), R[[j]]))
        
        Z_weighted = Z * weight_mat
        
        # MF_fixed
        
        MF_fixed = rowSums(X * weight_mat)
        
        # MF_random
        
        MF_random = rowSums(X^2 * weight_mat^2)
        
        T_fixed_mf_null[i] = T_fixed_mf_null[i] + sum(rowSums(Z_weighted)^2 * (abs(MF_fixed) > qnorm((1 - threshold/2))))
        
        T_random_mf_null[i] = T_random_mf_null[i] + sum(rowSums(Z_weighted^2) * 
                                    (sapply(length(MF_random), 
                                            function(q) davies(MF_random[q], lambda=weight_mat[q,]^2)$Qq) < threshold))
        
      }
    } else {
      
      nstudy = length(weight)
      
      weight_mat = matrix(rep(weight, ncol(R)),byrow=T,nrow=ncol(R))
      
      X = t(mvrnorm(nstudy, rep(0, ncol(R)), R))
      Z = t(mvrnorm(nstudy, rep(0, ncol(R)), R))
      
      Z_weighted = Z * weight_mat
      
      # MF_fixed
      
      MF_fixed = rowSums(X * weight_mat)
      
      # MF_random
      
      MF_random = rowSums(X^2 * weight_mat^2)
      
      T_fixed_mf_null[i] = sum(rowSums(Z_weighted)^2 * (abs(MF_fixed) > qnorm((1 - threshold/2))))
      
      T_random_mf_null[i] = sum(rowSums(Z_weighted^2) * 
                                  (sapply(length(MF_random), 
                                          function(q) davies(MF_random[q], lambda=weight_mat[q,]^2)$Qq) < threshold))
      
    }
    
  }
  
  return(c(sum(T_random_mf_null>=T_random_mf,na.rm=T), sum(T_fixed_mf_null>=T_fixed_mf,na.rm=T)))
}


## Fisher's method aggregate p-values from two tests

Fisher.test = function(p) {
  Xsq = -2*sum(log(p))
  p.val = pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
  return(list(Xsq = Xsq, p.value = p.val))
}
