#' Predictive value of core entries
#'
#' @param pheno pheno file
#' @param trait.col trait col
#' @param geno.rel.all geno matrix
#' @param core.names core names
#'
#' @return
#' @export
#'
#' @examples
pred.value<-function(pheno, trait.col, geno.rel.all, core.names){

test<-core.names
valid<-rownames(pheno)[!match(rownames(pheno),core.names, nomatch=0)>0] #validation PIs

#mask validation of phenotypes
train<-as.data.frame(cbind(rownames(pheno),pheno[,trait.col]))
train[(match(train$V1, valid, nomatch=0)>0),2]<-NA

#kin.blup
ans1<-rrBLUP::kin.blup(train, K=geno.rel.all, geno="V1", pheno="V2")
cor<-cor(ans1$g[valid], pheno[valid,trait.col], use="complete.obs")
return(cor)
}
