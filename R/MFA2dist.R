#' MFA to distance matrix
#'
#' This function calculates a distance matrix from a multiple factor analysis output.
#'
#' @param mfa a "MFA" object
#' @param n.dim numeric value specifying the number of principal components to use to construct a distance matrix
#' @param method the distance measure to be used. Default is "euclidean". See "stats::dist" for other options.
#'
#' @return MFA2dist.df returns a data frame of a distance object
#' @export
#'
MFA2dist.df<-function(mfa, n.dim, method="euclidean"){

  # if(rownames(geno.dist)!=rownames(pheno.dist)) stop('genotypic and phenotypic data are not in the same order')
  # dim.geno<-dim(as.data.frame(geno.dist))[2]
  # dim.pheno<-dim(as.data.frame(pheno.dist))[2]
  #
  # comb<-cbind(as.data.frame(geno.dist), as.data.frame(pheno.dist))
  #
  # mfa<-FactoMineR::MFA(comb, c(dim.geno, dim.pheno))

  return(as.data.frame(as.matrix(stats::dist(mfa$ind$coord[,1:n.dim], method=method))))

}
