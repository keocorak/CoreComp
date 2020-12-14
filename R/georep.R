#' Plot geographic map of core entries
#'
#' Plots geographic passport data for accessions
#'
#' @param x data frame containing columns for accessions names, latitude (column name must include "lat") and longitude (column name must include "long"). Coordinates must be in decimal form.
#' @param core.names character vector of accessions in core
#' @param core.groups vector or NULL: entry subgroup (likely from clustering). If provided, only core entries will be plotted
#' @param region character vector: regions/countries to plot. Default is a world map.
#' @param size integer: size of points to plot
#' @param ... other parameters
#'
#' @return map of accessions plotted according to geographic data. Core entries are either highlighted in red and all other accessions are plotted in black or core entries are colored according to subgroup and other accessions are not plotted.
#' @export
#'

plot.map<-function(x, core.names, core.groups=NULL, region=".", size=2,  ...){

  map<- NULL
  mapWorld <- ggplot2::borders("world", regions=region, colour="gray50", fill="gray50")
  map <- ggplot2::ggplot() + mapWorld



  if(is.null(core.groups)){
    print(map+ggplot2::geom_point(data=x, size=size,
                         ggplot2::aes_string(x=grep("long", colnames(x)), y=grep("lat", colnames(x)),
                                    color=factor(match(x[,1],core.names, nomatch=0)>0)),
                        position=ggplot2::position_jitter(w=10,h=10))+
            ggplot2::scale_color_manual(values=c("black", "red"))+
            ggplot2::labs(colour="core entry T/F")+
            ggplot2::theme_bw()+ggplot2::ylim(c(-60,85)))
  } else {
    print(map+ggplot2::geom_point(data=x[match(x[,1],core.names, nomatch=0)>0,], size=size,
                         ggplot2::aes_string(x=grep("long", colnames(x)), y=grep("lat", colnames(x)),
                                    color=factor(core.groups)),
                         position=ggplot2::position_jitter(w=10,h=10))+
            ggplot2::labs(colour="core cluster")+
            ggplot2::theme_bw()+ggplot2::ylim(c(-60,85)))
  }


}


