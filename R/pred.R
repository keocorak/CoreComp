#' Predictive value of core entries
#'
#' @param pheno dataframe of phenotypic information file
#' @param trait.col trait column index
#' @param geno.rel.all genotypic  relationship matrix matrix
#' @param core.names character vector of core entries
#'
#' @return correlation between predicted and actual value for accessions outside of the core
#' @export
#'


pred.accuracy<-function(pheno, trait.col, geno.rel.all, core.names){

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
