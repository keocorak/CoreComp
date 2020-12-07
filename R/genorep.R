#' Average distance between each accession and the nearest entry (A.NE)
#'
#' Calculate distance between each accession and the nearest entry in the core and average over all accessions
#'
#' @param geno.dist.all matrix
#' @param core.names vector
#'
#' @return
#' @export
#'
#' @examples
ANE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.ane<-as.matrix(geno.dist.all[core.names, !colnames(geno.dist.all)  %in% core.names]) #compare core to all
  l<-apply(geno.dist.ane,2,min) #find min
  sum(l)/length(l) #find average
}


#' Average distance between each entry and nearest neighbor entry  (E.NE)
#'
#' Calculate distance between nearest entries in the core and average over all entries
#'
#' @param geno.dist.all matrix
#' @param core.names vector
#'
#' @return
#' @export
#'
#' @examples
ENE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.ene<-as.matrix(geno.dist.all[core.names, core.names])
  diag(geno.dist.ene)<-NA
  s<-apply(geno.dist.ene, 2, min, na.rm=T)
  sum(s)/length(s)
}


#' Average distances between entries
#'
#' Sum distances between all pairs of entries in the core and divide by number of comparisons
#'
#' @param geno.dist.all matrix
#' @param core.names vector
#'
#' @return
#' @export
#'
#' @examples
EE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.ee<-as.matrix(geno.dist.all[core.names, core.names])
  geno.dist.ee[upper.tri(geno.dist.ee)]<-NA
  ee<-apply(geno.dist.ee, 2, sum, na.rm=T)
  #divide sum by number of pairwise compariosons
  sum(ee)/(dim(geno.dist.ee)[1]*(dim(geno.dist.ee)[1]-1)/2)
}
