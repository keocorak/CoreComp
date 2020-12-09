
#'Correlate Cophenetic Distances
#'
#' calculate correlation between cophenetic distances and original distance matrix
#'
#' @param geno.dist.all genotypic distance matrix
#' @param hc_func agglomeration method to be used. Default method is "ward.D2" See stats::hclust for other options.
#'
#' @return
#' @export
#'
cor.cophentic<-function(geno.dist.all, hc_func="ward.D2"){
  cor(as.dist(geno.dist.all), cophenetic(hclust(as.dist(geno.dist.all), method=hc_func)))
}


sil.width<-function(geno.dist.all, n.clust=9, hc_func="hclust", hc_method="ward.D2", plot.sil=TRUE){


  if(n.clust>=dim(geno.dist.all)[1]) stop('number of clusters must be less than number of entries')


  sil_width<-as.data.frame(matrix(NA, nrow=n.clust-1, ncol=2))
  colnames(sil_width)<-c("n_cluster", "ave_sil_width")

    sil_width[,1]<-2:n.clust

  for (j in 2:n.clust){
    tmp<-factoextra::hcut(as.dist(geno.dist.all), k=j, isdiss=TRUE, hc_func=hc_func, hc_method=hc_method)
    sil_width[(j-1),2]<-tmp$silinfo$avg.width
  }


  sil.plot<- ggplot2::ggplot(data=sil_width, ggplot2::aes(x=n_cluster, y=ave_sil_width))+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept=sil_width[which.max(sil_width[,2]),1], linetype="dotted")

  if(plot.sil==TRUE){
    print(sil.plot)
  }

  return(sil_width)

}

# geno.core<-function(){
#
# }
