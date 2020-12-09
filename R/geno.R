
#'Correlate Cophenetic Distances
#'
#' calculate correlation between cophenetic distances and original distance matrix
#'
#' @param geno.dist.all genotypic distance matrix
#' @param hc_method agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#'
#' @return
#' @export
#'
cor.cophentic<-function(geno.dist.all, hc_method="ward.D2"){
  cor(as.dist(geno.dist.all), cophenetic(hclust(as.dist(geno.dist.all), method=hc_func)))
}


#' Calculate Average Silhouette Widths
#'
#' @param geno.dist.all genetic distance matrix
#' @param n.clust number of clusters to test, must be less than number of accessions
#' @param hc_func hierarchial clustering function to be used. Default is "hclust". See factoextra::hcut for other options.
#' @param hc_method agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#' @param plot.sil TRUE/FALSE to plot silhouette widths
#'
#' @return
#' @export
#'
#' @examples
sil.width<-function(geno.dist.all, n.clust.to.check=9, hc.func="hclust", hc.method="ward.D2", plot.sil=TRUE){


  if(n.clust>=dim(geno.dist.all)[1]) stop('number of clusters must be less than number of accessions')


  silwidth<-as.data.frame(matrix(NA, nrow=n.clust-1, ncol=2))
  colnames(silwidth)<-c("n.clusters", "ave.sil.width")

    silwidth[,1]<-2:n.clust.to.check

  for (j in 2:n.clust){
    tmp<-factoextra::hcut(as.dist(geno.dist.all), k=j, isdiss=TRUE, hc.func=hc.func, hc.method=hc.method)
    sil.width[(j-1),2]<-tmp$silinfo$avg.width
  }


  sil.plot<- ggplot2::ggplot(data=sil.width, ggplot2::aes(x=n.clusters, y=ave.sil.width))+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept=sil.width[which.max(sil.width[,2]),1], linetype="dotted")

  if(plot.sil==TRUE){
    print(sil.plot)
  }

  return(sil.width)

}

geno.core<-function(geno.dist.all, n.clust, hc.method="ward.D2", size.core=0.1){

    geno.clust<-hclust(as.dist(geno.dist.all), method=hc.method)
    group<-cutree(geno.clust, k=n.clust)


    group<-cbind(read.table(text=names(group)), group)
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


