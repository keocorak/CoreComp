
#'Correlate Cophenetic Distances
#'
#' Calculates correlation between cophenetic distances and an original distance matrix.
#' High correlation (>0.8) can provide evidence for population structure.
#'
#' @param dist.all distance matrix
#' @param hc.method agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#'
#' @return correlation between cophenetic distances and original distance matrix
#' @export
#'
#' @examples data(fake_geno); cor.cophenetic(dist(fake_geno))
cor.cophenetic<-function(dist.all, hc.method="ward.D2"){
  stats::cor(stats::as.dist(dist.all), stats::cophenetic(stats::hclust(stats::as.dist(dist.all), method=hc.method)))
}


#' Calculate Average Silhouette Widths
#'
#' Returns and plots average silhouette widths for 1:n clusters where n is any number less than the total number of accessions.
#'
#' @param dist.all distance matrix
#' @param n.clust.to.check integer function will evaluate cluster sizes from 1:n.clust.to.check
#' @param hc.func hierarchical clustering function to be used. Default is "hclust". See factoextra::hcut for other options.
#' @param hc.method agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#' @param plot.sil logical TRUE/FALSE to output scree plot of silhouette widths
#'
#' @return data frame of average silhouette width for different numbers of clusters
#' @export
#'
#' @examples
#' data(fake_geno); sil.width(dist(fake_geno))

sil.width<-function(dist.all, n.clust.to.check=9, hc.func="hclust", hc.method="ward.D2", plot.sil=TRUE){


  if(n.clust.to.check>=dim(as.matrix(dist.all))[1]) stop('number of clusters must be less than number of accessions')


  silwidth<-as.data.frame(matrix(NA, nrow=n.clust.to.check-1, ncol=2))
  colnames(silwidth)<-c("n.clusters", "ave.sil.width")

    silwidth[,1]<-2:n.clust.to.check

  for (j in 2:n.clust.to.check){
    tmp<-factoextra::hcut(stats::as.dist(dist.all), k=j, isdiss=TRUE, hc.func=hc.func, hc.method=hc.method)
    silwidth[(j-1),2]<-tmp$silinfo$avg.width
  }


  sil.plot<- ggplot2::ggplot(data=silwidth, ggplot2::aes_string(x='n.clusters', y='ave.sil.width'))+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept=silwidth[which.max(silwidth[,2]),1], linetype="dotted")

  if(plot.sil==TRUE){
    print(sil.plot)
  }

  return(silwidth)

}

#' Select a Core Collection based on a distance matrix
#'
#' Clusters accessions into groups based on pairwise distances and takes a random subset from each cluster to compose a core set.
#'
#' @param dist.all a distance matrix
#' @param n.clust number of clusters into which to group accessions
#' @param size.core percentage of total number of accessions to select for core, default is 10\%
#' @param hc.method agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#'
#' @return dataframe of core entries and the cluster group, pass entry column as a vector
#' @export
#'
#' @examples
#' data(fake_geno); dist.core(dist(fake_geno), n.clust=9)


dist.core<-function(dist.all, n.clust, size.core=0.1, hc.method="ward.D2"){

  if(n.clust>=dim(as.matrix(dist.all))[1]) stop('number of clusters must be less than number of accessions')
  if(size.core>1) stop('size of core must be expressed as a proportion')

    dist.clust<-stats::hclust(stats::as.dist(dist.all), method=hc.method)
    group<-stats::cutree(dist.clust, k=n.clust)


    group<-cbind(utils::read.table(text=names(group)), group)
    colnames(group)[1]<-"name"


    group.list<-split(group, group$group)

    group.list.core<-vector("list", n.clust)

    for(i in 1:n.clust){
      x<-group.list[[i]][sample(nrow(group.list[[i]]),round(size.core*dim(group.list[[i]])[1])),]
      group.list.core[[i]]<-x
    }


    group.list.core<-plyr::ldply(group.list.core, data.frame)
    group.list.core$group<-factor(group.list.core$group, levels=c(1:n.clust))

    return(group.list.core)

  }


