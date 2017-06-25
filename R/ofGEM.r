#' A meta-analysis approach with filtering for identifying gene-level gene-environment interactions with genetic association data
#' 
#'This function first conducts a meta-filtering test to filter out unpromising SNPs. It then runs a test of omnibus-filtering-based GxE meta-analysis (ofGEM) that combines the strengths of the fixed- and random-effects meta-analysis with meta-filtering. It can also analyze data from multiple ethnic groups. The p-values are calculated using a sequential sampling approach.
#' 
#' 
#' @param Z the test statistics for gene-environment interactions (GxE). Each row is a SNP, and each column is a study. For multi-ethnic groups,
#' Z is a list with each element as a matrix for each ethnic group.
#' @param X the filtering statistics for GxE. Each row is a SNP, and each column is a study. For multi-ethnic groups,
#' X is a list with each element as a matrix for each ethnic group.
#' @param R the correlation matrix to simulate test and filtering statistics under the null distribution. 
#' The simulated test and filtering statistics are used for testing. For multi-ethnic groups,
#' R is a list with each element as a correlation matrix for each ethnic group.
#' @param weight the weight vector for each study, or the weight matrix for each SNP and each study. If the weight is common
#' across SNPs, it is a vector with a length equal to the number of studies. If the weight is different across SNPs, it is a matrix
#' with each row for a SNP and each column as a study.
#' @param threshold a fixed p-value threshold for filtering. The default is 0.1.
#' @param maxSim the maximum number of samples to be simulated for the test and filtering statistics under the null. The default is 1e6.
#' @param tol the tolerance number to stop the sequential sampling procedure. We count the number of simulated test statistics
#' with an absolute value larger than that of the calculated one based on the data, for every 100 simulations. 
#' The sampling will stop if the count reaches tol. The default is 10. If tol = 0, the number
#' of simulations equals to maxSim.
#' 
#' @return A list containing
#' \item{pval_random_mf}{the p-value based on random-effects meta-analysis with meta-filtering.}
#' \item{pval_fixed_mf}{the p-value based on fixed-effects meta-analysis with meta-filtering.}
#' \item{pval_ofGEM}{the p-value based on aggregating the p-values of fixed- and random-effects meta-analyses with meta-filtering 
#' using Fisher's method.}
#' \item{nsim}{the number of simulations that are performed.}
#' 
#' @references Wang, Liu, Pierce, Huo, Nicolae, Olopade, Ahsan, & Chen (2017+). A meta-analysis approach with filtering for 
#' identifying gene-level gene-environment interactions with genetic association data. In preparation.
#' 
#' @export
#' @import MASS CompQuadForm
#' @importFrom stats pchisq qnorm
#' 
#' @examples
#' 
#' data(sim_dat)
#' 
#' pval = ofGEM(Z = sim_dat$Z, X = sim_dat$X, R = sim_dat$R, weight = rep(1/6, 6))
#' 


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
    
    while (!all(count >= 10, na.rm = T) & nsim < maxSim) {
      
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
