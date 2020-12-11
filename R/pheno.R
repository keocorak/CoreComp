#' Least Squares Means
#'
#' Calculate least squares means for experimental phenotype data
#'
#' @param x data table
#' @param id.vars numeric vector corresponding to columns in data
#' @param measure.vars numeric vector corresponding to columns in data
#' @param form model formula, for example "value~-1+entryID+(1|rep)"
#' @param uniquecol integer corresponding to entry names column- will be included in data frame output
#' @param like character string that indicates that something is an entry
#' @param alpha confidence level for confidence intervals
#'
#' @return Rdata files of lsmeans and confidence intervals
#' @export
#'
#'
lsmeans<-function(x,id.vars, measure.vars, form, uniquecol, like,  alpha=0.05){

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
    blues<-lme4::fixef(fit, add.dropped=T)[grep(like, names(lme4::fixef(fit)))]
    blues.frame[levels(fit@frame[[2]]), trait]<-blues

    se_blues <- stats::coef(summary(fit))[grep(like, names(lme4::fixef(fit))), "Std. Error"]
    edf <- sum(stats::hatvalues(fit))
    ci <- matrix(rep(blues, 2), ncol=2) + (stats::qt(p=1-alpha/2, df=edf) * se_blues) %*% t(c(-1, 1))
    ci.array[levels(fit@frame[[2]]), trait, ] <- ci
  }
  saveRDS(blues.frame, file="lsmeans.Rdata")
  saveRDS(ci.array, file="ci.Rdata")
}
