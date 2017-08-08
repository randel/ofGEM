

#' A simulated data example
#' 
#' This simulated data list is for demonstration.
#' 
#' 
#' @name sim_dat
#' @docType data
#' @return A list containing \item{Z}{a simulated matrix of testing statistics
#' for GxE. Each row corresponds to a SNP in a set, and each column represents
#' a study.} \item{X}{a simulated matrix of filtering statistics for GxE. Each
#' row corresponds to a SNP, and each column represents a study.} \item{R}{a
#' simulated correlation matrix for testing (or filtering) statistics in a set
#' under the null. The testing statistics for SNPs in LD will be also
#' correlated. The correlations can be approximated by the LD matrix. We
#' resample SNP-level testing and filtering statistics from multivariate normal
#' distributions to calculate the null meta-analysis statistics.}
#' @examples
#' 
#' data(sim_dat)
#' 
NULL



