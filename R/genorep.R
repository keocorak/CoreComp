#' Average distance between each accession and the nearest entry (A.NE)
#'
#' Calculate distance between each accession and the nearest entry in the core and then average over all accessions. Value should be as small as possible (see Odong et al. 2013).
#'
#' @param geno.dist.all genetic distance matrix
#' @param core.names character vector of entry names to be included in core
#'
#' @return average distance between each accession and the nearest entry
#' @export
#'
#' @examples core<-dist.core(as.matrix(dist(fake_geno)), n.clust=9); ANE(dist(fake_geno), core$name)
#'
ANE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.all<-as.matrix(geno.dist.all)
  geno.dist.ane<-as.matrix(geno.dist.all[as.vector(core.names), !colnames(geno.dist.all)  %in% as.vector(core.names)]) #compare core to all
  l<-apply(geno.dist.ane,2,min) #find min
  sum(l)/length(l) #find average
}


#' Average distance between each entry and nearest neighbor entry  (E.NE)
#'
#' Calculate distance between nearest entries in the core and average over all entries. Value should be as large as possible (see Odong et al. 2013).
#'
#' @param geno.dist.all genetic distance matrix
#' @param core.names character vector of entry names to be included in core
#'
#' @return average distance between each entry and nearest neighbor entry
#' @export
#'
#' @examples core<-dist.core(as.matrix(dist(fake_geno)), n.clust=9); ENE(dist(fake_geno), core$name)
ENE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.all<-as.matrix(geno.dist.all)
  geno.dist.ene<-as.matrix(geno.dist.all[as.vector(core.names), as.vector(core.names)])
  diag(geno.dist.ene)<-NA
  s<-apply(geno.dist.ene, 2, min, na.rm=T)
  sum(s)/length(s)
}


#' Average distances between entries (EE)
#'
#' Sum distances between all pairs of entries in the core and divide by number of comparisons. Value should be as large as possible (see Odong et al. 2013).
#'
#' @param geno.dist.all distance matrix
#' @param core.names character vector of entry names to be included in core
#'
#' @return average distance between entries
#' @export
#'
#' @examples core<-dist.core(as.matrix(dist(fake_geno)), n.clust=9); ENE(dist(fake_geno), core$name)
EE<-function(geno.dist.all, core.names){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  geno.dist.all<-as.matrix(geno.dist.all)
  geno.dist.ee<-as.matrix(geno.dist.all[as.vector(core.names), as.vector(core.names)])
  geno.dist.ee[upper.tri(geno.dist.ee)]<-NA
  ee<-apply(geno.dist.ee, 2, sum, na.rm=T)
  #divide sum by number of pairwise comparisons
  sum(ee)/(dim(geno.dist.ee)[1]*(dim(geno.dist.ee)[1]-1)/2)
}


#' Noirot Principal Component scoring
#'
#' Perform multidimensional scaling of genotypic distance matrix and calculate summed relative contribution of entries in core set following Noirot et al. (1996).
#'
#' @param geno.dist.all genetic distance matrix
#' @param core.names character vector of entry names to be included in core
#' @param k the maximum dimensions in which to represent data; must be less than n
#'
#' @return PC score
#' @export
#'
#' @examples core<-dist.core(as.matrix(dist(fake_geno)), n.clust=9)
#' @examples noirot.contribution(as.matrix(dist(fake_geno)), core$name)

noirot.contribution<-function(geno.dist.all, core.names, k=2){
  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  if(k>dim(geno.dist.all)[1]) stop('k must be less than n')

  geno.dist.all.mds<-stats::cmdscale(geno.dist.all,k)

  cum.all<-list()
  for(j in 1:dim(geno.dist.all.mds)[1]){
    all<-list()
    for(i in 1:dim(geno.dist.all.mds)[2]){
      all[i]<-(geno.dist.all.mds[j,i])^2
    }
    cum.all[j]<-sum(unlist(all))
  }

  total<-sum(unlist(cum.all))



  geno.dist.all.mds.core<-geno.dist.all.mds[rownames(geno.dist.all.mds)%in% as.vector(core.names), ]

  jsum<-list()
  for(j in 1:dim(geno.dist.all.mds.core)[1]){ #for every entry in k

    rowsum<-list()
    for(i in 1:dim(geno.dist.all.mds.core)[2]){ #square each coordinate
      rowsum[i]<-(geno.dist.all.mds.core[j,i]^2)
    }

    jsum[j]<-sum(unlist(rowsum)) #and sum it

  }

  sum(unlist(jsum))/total


  }


#' MDS bi-plot of accessions
#'
#' @param x genetic distance matrix
#' @param core.names character vector of entry names to be included in core
#' @param k the maximum dimensions in which to represent data; must be less than n
#' @param axes vector of length 2 specifying dimensions to represent with bi-plot
#' @param ... other parameters
#'
#' @return biplot of results from multidimensional scaling of genetic distance matrix with entries selected for the core highlighted
#' @export
#'
plot.geno.entries<-function(x, core.names, k=2, axes=c(1,2), ...){

  geno.dist.all=x

  if(is.null(core.names)) stop('provide vector of entries')
  if(is.null(geno.dist.all)|is.vector(geno.dist.all)) stop('matrix cannot be a vector')
  if(k>dim(geno.dist.all)[1]) stop('k must be less than n')

  geno.dist.all.mds<-stats::cmdscale(geno.dist.all,k)
  ggplot2::ggplot(as.data.frame(geno.dist.all.mds), ggplot2::aes(x=geno.dist.all.mds[,axes[1]], y=geno.dist.all.mds[,axes[2]]))+
    ggplot2::geom_point(ggplot2::aes(color=factor(rownames(geno.dist.all.mds) %in% as.vector(core.names))))+
    ggplot2::xlab(paste("PC", axes[1], sep=" "))+ggplot2::ylab(paste("PC", axes[2], sep=" "))+
    ggplot2::scale_color_discrete(name="core entries", labels=c("no", "yes"))


  }








