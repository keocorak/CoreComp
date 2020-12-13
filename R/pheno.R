#' Least Squares Means
#'
#' Calculate least squares means for experimental phenotype data
#'
#'
#'
#' @param x a data frame containing entry names,
#' @param id.vars numeric vector corresponding to ID variables columns in data
#' @param measure.vars numeric vector corresponding to measured variables columns in data
#' @param form model formula, for example "value~-1+entryID+(1|rep)"
#' @param uniquecol integer corresponding to entry names column- will be included in data frame output
#' @param out data to return. Should be either "means" for least-squares means estimates or "ci" for confidence intervals.
#' @param alpha significance level used to construct confidence intervals. Only used when out="ci"
#' @return dataframe of means or confidence intervals
#' @export
#'
#'
lsmeans<-function(x,id.vars, measure.vars, form, uniquecol, out="means", alpha=NULL){

  #thanks to guilliame ramstein for most of this code. I just put it into a function

  meltdata<-reshape::melt(x, id.vars=names(x[id.vars]), measure.vars=names(x[measure.vars]), value.name="value")
  #meltdata[id.vars]<-lapply(meltdata[id.vars],factor)
  traits<-unique(meltdata$variable)
  blues.frame<-matrix(numeric(),
                      nrow=length(unique(meltdata[,uniquecol])), ncol=length(unique(meltdata$variable)),
                      dimnames=list(unique(meltdata[,uniquecol]), unique(meltdata$variable)))
  ci.array<- array(numeric(),
                   dim=c(length(unique(meltdata[,uniquecol])), length(unique(meltdata$variable)), 2),
                   dimnames=list(unique(meltdata[,uniquecol]), unique(meltdata$variable), c("L","U"))
  )

  for(trait in traits){
    tmp<-meltdata[meltdata$variable==trait,]
    tmp<-tmp[!is.na(tmp$value),]
    tmp$value<-as.numeric(tmp$value)
    form<-stats::as.formula(form)
    fit<-lme4::lmer(form, data=tmp)
    blues<-lme4::fixef(fit, add.dropped=T)[grep(colnames(x)[uniquecol], names(lme4::fixef(fit)))]
    blues.frame[levels(fit@frame[[2]]), trait]<-blues
  }

  if(out=="means"){
    return(blues.frame)
  } else if(out=="ci") {
    se_blues <- stats::coef(summary(fit))[grep(colnames(x)[uniquecol], names(lme4::fixef(fit))), "Std. Error"]
    edf <- sum(stats::hatvalues(fit))
    ci <- matrix(rep(blues, 2), ncol=2) + (stats::qt(p=1-alpha/2, df=edf) * se_blues) %*% t(c(-1, 1))
    ci.array[levels(fit@frame[[2]]), trait, ] <- ci
    return(ci.array)
  }
}




