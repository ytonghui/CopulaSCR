#################################################################
##   CopulaSCR R package by Tonghui Yu Copyright (C) 2024
##
## The below functions are used for survival prediction for Semi-competing
## risks data with single intermediate event times
#################################################################
#' @title Survival prediction
#' @description
#' Estimation of  joint survival function, cross-hazard ratio function
#' and conditional survival function under semi-competing risks data with
#' single intermediate event.
#'
#' @name  predictscr
#' @aliases predictscr
#' @aliases predjsurv
#' @aliases predcsurv
#' @aliases predchr
#'
#' @usage predictscr(data, fit, t1.formula, t2.formula, type = c("jsurv","csurv","chr"),
#'    t1, t2, s2, t1equal = FALSE, method = c("np", "sp"),...)
#' @usage predjsurv(data, fit, t1.formula, t2.formula, t1, t2, method = c("np", "sp"),newdata)
#' @usage predchr(data,fit,t1.formula, t2.formula, t1, t2, method=c("np","sp"),
#'    tau, copulafam = c("clayton","frank","joe","gumbel","amh"),newdata)
#' @usage predcsurv(data, fit, t1.formula, t2.formula, t1, t2, s2,
#'   t1equal = FALSE, method = c("np","sp"), newdata)
#'
#' @param data an optional data frame in which to interpret the variables
#'   occurring in the formula.
#' @param fit an optional object of class  "scrsurv".
#' @param type specifies the type of function to be estimated.
#' @param t1.formula  an optional formula expression for the nonterminal event
#'   time, of the form \code{response ~ predictors}. The \code{response} is a
#'   \code{Surv} object with right censoring. See the documentation of
#'   \code{coxph} and \code{formula} for details.
#' @param t2.formula  an optional formula expression for the terminal event
#'   time, of the form \code{response ~ predictors}.
#' @param t1 numeric vector
#' @param t2 numeric vector
#' @param s2 numeric vector
#' @param t1equal If TRUE,..
#' @param method a character string specifying the method for computing survival
#'   curves.
#' @param tau An optional numeric value of Kendall's tau parameter in the Archimedean copula.
#' @param copulafam An optional character string specifying the family of an Archimedean
#'   copula.
#' @param newdata input new data frame
#' @param ... other arguments
#'
#' @return
#' \describe{
#' \item{jsurv}{A matrix for joint survival function \eqn{Pr(\tilde{T}_1 > t_1, \tilde{T}_2 > t_2)},
#'  where \eqn{\tilde{T}_1}, and \eqn{\tilde{T}_2}
#'    are the true nonterminal and terminal event times, respectively.}
#' \item{chr}{A matrix for cross-hazard ratio function.}
#' \item{csurv}{A matrix for conditional survival function
#' \eqn{Pr(\tilde{T}_2 > t_2|\tilde{T}_1 > t_1, \tilde{T}_2 > s_2)} if  t1equal = FALSE
#' and \eqn{Pr(\tilde{T}_2 > t_2|\tilde{T}_1 = t_1, \tilde{T}_2 > s_2)} if t1equal = TRUE. }
#' }
#' \item{...}{}
#' @export predictscr
#' @export predjsurv
#' @export predchr
#' @export predcsurv
#'
#' @seealso \code{\link{scrassonp}}, \code{\link{scrsurv}},\code{\link{predict.scrsurv}}
#' @examples library(CopulaSCR)
#' set.seed(12345)
#' simdata<- simSCR(n = 50,tau = 0.5, copulafam = "frank")
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
#'                     t2.formula=Surv(T2, event2) ~ 1,
#'                     data = simdata,copulafam="frank",
#'                     a=quantile(simdata$T1,0.9),
#'                     b=quantile(simdata$T2,0.9),B=10,seed = 12345,
#'                     se = TRUE,se.method = "resampling")
#' summary(fitasso)
#' fit22<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=fitasso$B,
#' se.method = fitasso$se.method,conf.int=TRUE,conftype =3)
#'
#' t1<- t2<- seq(0,3.5,0.1)
#' ### joint survival probabilities
#' js<- predjsurv(fit=fit22,t1=t1, t2=t2,method="sp")
#' # same with
#' js2<- predictscr(fit=fit22,t1=t1, t2=t2,method="sp",type="jsurv")
#' # or semiparametric estimator
#' # js2<- predict(fit=fit22,t1=t1, t2=t2,type="jsurv")
#'
#' ### cross-hazard ratio function
#' chr<- predchr(fit=fit22, t1=t1, t2=t2,method="sp")
#' # same with
#' chr2<- predictscr(fit=fit22,t1=t1, t2=t2,method="sp",type="chr")
#' # or semiparametric estimator
#' # chr2<- predict(fit=fit22,t1=t1, t2=t2,type="chr")
#'
#'
#' ## plot jsurv and chr
#' library(graphics)
#' filled.contour(x=t1,y=t2,z=do.call(rbind,js$jsurv),nlevels = 30,
#' xlab="t1",ylab="t2",main = "joint survival function",
#' plot.axes = {axis(1, seq(0,3.5,0.5))
#'   axis(2, seq(0,3.5,0.5))
#'     contour(x = t1,y=t2,do.call(rbind,js$jsurv),
#'     add = TRUE, lwd = 2)})
#'
#'
#' filled.contour(x=t1,y=t2,z=do.call(rbind,chr$chr),nlevels = 30,
#' xlab="t1",ylab="t2",main = "cross-hazard ratio",
#' plot.axes = {axis(1, seq(0,3.5,0.5))
#'  axis(2, seq(0,3.5,0.5))
#'  contour(x = t1,y=t2,do.call(rbind,chr$chr),
#'  add = TRUE, lwd = 2)})
#'
#' # conditional survival probabilities
#' t1<- s2<- c(0.2,0.5,1)
#' t2<-  seq(0,3.5,0.1)
#' cs<- predcsurv(fit=fit22,t1=t1, t2=t2, s2=s2, t1equal =TRUE,method="sp")$csurv
#' # same with
#' cs2<- predictscr(fit=fit22,t1=t1, t2=t2, s2=s2, t1equal =TRUE,
#'  method="sp",type="csurv")$csurv
#' # or semiparametric estimator
#' # cs2<- predict(fit=fit22,t1=t1, t2=t2,s2=s2, t1equal =TRUE,type="csurv")$csurv
#'
#' ## plot csurv
#' cs<- lapply(seq(length(cs)),function(k)cs[[k]][,k])
#' plot(c(0,max(t2)),c(0,1),type="n",xlab = "t2", ylab="",
#'   main ="conditional survival given T1=t1, T2>t1")
#' for (k in 1:length(cs)) {lines(t2,cs[[k]],col=k)}
#' legend(x="topright",legend = paste0("t1=",t1),col = 1:length(cs),lty=1)
#'
#'



predictscr<- function(data, fit, t1.formula, t2.formula, type = c("jsurv","csurv","chr"),
                      t1, t2, s2, t1equal = FALSE, method = c("np", "sp"), ...){
  Call <- match.call()
  method<- match.arg(method)
  type <- match.arg(type)
  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2' ,'s2','tau','copulafam'),
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  if (type=="csurv"&indx[7]==0) stop("a s2 argument is required")


  if(method=="np"){
    if (indx[1]==0&indx[2]==0) stop("a data argument is required")
    if (indx[3]==0&indx[2]==0) stop("a t1.formula argument is required")
    if (indx[4]==0&indx[2]==0) stop("a t2.formula argument is required")
  }


  if(indx[1]==0) data<- NULL
  if(indx[2]==0) fit<- NULL
  if(indx[3]==0) t1.formula<- NULL
  if(indx[4]==0) t2.formula<- NULL

  if(type=="jsurv"){
    out<- predjsurv(data=data, fit=fit, t1.formula =t1.formula, t2.formula =t2.formula,
                    t1=t1, t2=t2, method = method, ...)
  }
  if(type=="csurv"){
    out<- predcsurv(data=data, fit=fit, t1.formula =t1.formula, t2.formula =t2.formula,
                    t1=t1, t2=t2, s2=s2, t1equal = t1equal, method = method, ...)
  }
  if(type=="chr"){
    out<- predchr(data=data, fit=fit, t1.formula =t1.formula, t2.formula =t2.formula,
                  t1=t1, t2=t2, method = method,...)
  }
  return(out)
}

predjsurv<- function(data,fit,t1.formula, t2.formula,t1, t2,method=c("np","sp"),newdata){
  Call <- match.call()
  method<- match.arg(method)

  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2'),
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  # if(length(t1)!=length(t2)) stop("t1 and t2 shall have equal length")

  if(method=="np"){
    if(indx[2]==0){
      if (indx[1]==0) stop("a data argument is required")
      if (indx[3]==0) stop("a t1.formula argument is required")
      if (indx[4]==0) stop("a t2.formula argument is required")
    }


    if(indx[2]!=0) {
      if(is.null(fit)){stop("fit shall be 'scrsurv' class.")}
      if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}
      t1.formula<- fit$t1.formula
      t2.formula<- fit$t2.formula
      if(indx[1]!=0){if(is.null(data))data<- fit$data}
    }
    if(indx[1]!=0){
      data<- prepData0(t1.formula,t2.formula, data)
      na.action<- data$na.action
      data<- data$data
      t1data <-  prepData(formula=t1.formula, data, na.action,bandwidth, discrete.level,xlimits)
      t2data <-  prepData(formula=t2.formula, data, na.action,bandwidth, discrete.level,xlimits)


    }else{
      t1data<- fit$t1data
      t2data<- fit$t2data

    }

    if(length(all.vars(t1.formula[[3]]))>1)
      {stop("nonparametric estimation for stratified joint survival function is under construction.")}
    time1<- t1data$time
    event1<- t1data$event
    time2<- t2data$time
    event2<- t2data$event
    n<- length(time1)

    js0<- lapply(1:length(t1),function(k1)
      sapply(1:length(t2),function(k2)
        mean((time1>t1[k1])*(time2>t2[k2]))
      ))


    if(sum(1-event2)/n>0.05){
      Gfit<- survival::survfit(survival::Surv(time2,1-event2)~1)
      Gest.t2<- approx(c(0,Gfit$time),c(1,Gfit$surv),xout =  t2,method="constant",f = 0)$y
      rm(Gfit)

      js<- lapply(js0,function(xx)ifelse(xx==0& Gest.t2==0,0, xx/ Gest.t2))

    }else{
      js<- js0
    }
    names(js)<-  paste0("t1=",t1)
    return(list(jsurv=js,t1=t1,t2=t2))
  }

  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if(is.null(fit)){stop("fit shall be 'scrsurv' class.")}
    if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}

    surv1<- fit$t1.surv
    surv2<- fit$t2.surv
    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam

    if (fit$t1.surv$covariate.type==1){
      s1.t1<- approx(x=c(0,surv1$time),y=c(1,surv1$surv),xout = t1,method = "constant",
                     yleft = 1,yright = min(surv1$surv),f = 0,na.rm = TRUE)$y

      s2.t2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = t2,method = "constant",
                     yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE)$y
      js<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.t2))
      names(js)<-  paste0("t1=",t1)

    }else{
      strata.vars <- fit$t1.surv$discrete.predictors
      X.formula <- update(formula(fit$t1.formula),NULL~.)
      if (!all(match(all.vars(X.formula),names(newdata),nomatch=FALSE)))
        stop("Arg newdata does not contain all the covariates used for fitting. \n\nfitted variables: ", paste(all.vars(X.formula),collapse=", "),"\nnewdata contains:",ifelse(length(names(newdata))==0," nothing",names(newdata)))
      requested.X <- newdata[,all.vars(X.formula),drop=FALSE]


      fit.strata <- interaction(fit$stratas[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      requested.strata <- interaction(requested.X[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      fit.levels <- as.character(unique(fit.strata))
      NS <- length(fit.levels)
      fit.strata <- factor(fit.strata,levels=levels(fit.strata),labels=1:NS)
      requested.strata <- factor(requested.strata,levels=fit.levels,labels=1:NS)
      strata<- rep(fit.strata,times=fit$t1.surv$size.strata)
      copulaparam<- copulaparam[requested.strata]


      survmat1<- data.frame(time=surv1$time,surv=surv1$surv,strata=strata)
      survfn1<- lapply(fit.strata,function(k){
        sv<- survmat1[survmat1[,3]==k,]
        stats::stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
      })
      psurv1<- survfn1[requested.strata]
      s1.t1<- lapply(psurv1,function(x)x(t1))


      survmat2<- data.frame(time=surv2$time,surv=surv2$surv,strata=strata)
      survfn2<- lapply(fit.strata,function(k){
        sv<- survmat2[survmat2[,3]==k,]
        stats::stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
      })
      psurv2<- survfn2[requested.strata]
      s2.t2<- lapply(psurv2,function(x)x(t2))


      js<- lapply(1:length(s1.t1),function(k)
        {a<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param =copulaparam[k] ,
                 p1 = s1.t1[[k]][k1],p2 = s2.t2[[k]]))
        names(a)<-  paste0("t1=",t1)
        a
        })

      names(js) <- paste(strata.vars,requested.X[,1],sep="=")

    }


    return(list(jsurv=js,t1=t1,t2=t2,copulafam=copulafam,copulaparam=copulaparam))
  }
}


predchr<- function(data,fit,t1.formula, t2.formula, t1, t2,method=c("np","sp"),
                   tau,copulafam=c("clayton","frank","joe","gumbel","amh"),newdata){

  Call <- match.call()
  method<- match.arg(method)
  copulafam<- match.arg(copulafam)

  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2', 'tau'),
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")


  if(method=="np"){
    if(indx[2]==0){
      if (indx[1]==0) stop("a data argument is required")
      if (indx[3]==0) stop("a t1.formula argument is required")
      if (indx[4]==0) stop("a t2.formula argument is required")
      if (indx[7]==0) stop("a tau argument is required")
    }

    if(indx[2]!=0) {
      if(is.null(fit)) stop("fit shall be 'scrsurv' class.")

      if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}


      js<- predjsurv(fit = fit, t1=t1, t2=t2,method="np")
      copulafam<- fit$copulafam
      copulaparam<- fit$copulaparam

    }else{
      if (is.null(copulafam)) stop("a copulafam argument is required")
      if (is.null(tau)) stop("a tau argument is required")
      if (is.null(data)) stop("a data argument is required")
      if (is.null(t1.formula)) stop("a t1.formula argument is required")
      if (is.null(t2.formula)) stop("a t2.formula argument is required")

      copulaparam<- Calitau(copulafam = copulafam, tau = tau)
      js<- predjsurv(data = data,t1.formula = t1.formula, t2.formula = t2.formula,
                     t1=t1, t2=t2,method="np")

      chr<- lapply(js$jsurv,function(xx)
        npCHRfit(s=xx,param = copulaparam,copulafam = copulafam))
    }

    return(list(chr=chr,t1=t1,t2=t2))
  }

  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if (is.null(fit)) if (indx[2]==0) stop("a fit argument is required")
    if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}


    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam


    if (fit$t1.surv$covariate.type==1){
      js<- predjsurv(fit = fit, t1=t1, t2=t2,method="sp")
      chr<- lapply(js$jsurv,function(xx)
        npCHRfit(s=xx,param = copulaparam,copulafam = copulafam))

    }  else{
      js<- predjsurv(fit = fit, t1=t1, t2=t2,method="sp",newdata=newdata)
      copulaparam<- js$copulaparam
      js<- js$jsurv
      chr<- lapply(1:length(js),function(xx)
        npCHRfit(s=do.call(cbind,js[[xx]]),
                 param = copulaparam[xx],copulafam = copulafam))
      chr<- lapply(1:length(chr),function(xx){
        do.call(c,apply(chr[[xx]],2,list))
      })
      names(chr)<- names(js)
    }

    return(list(chr=chr,t1=t1,t2=t2,copulaparam=copulaparam,copulafam = copulafam))
  }






}



predcsurv<- function(data,fit,t1.formula, t2.formula,t1, t2, s2,
                       t1equal =FALSE,method=c("np","sp"),newdata){
  Call <- match.call()
  method<- match.arg(method)

  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2','s2'),
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  if (indx[7]==0) stop("a s2 argument is required")

  if(method=="np"){
    if(indx[2]==0){
      if (indx[1]==0) stop("a data argument is required")
      if (indx[3]==0) stop("a t1.formula argument is required")
      if (indx[4]==0) stop("a t2.formula argument is required")
    }


    if(indx[2]!=0) {
      if(is.null(fit)){stop("fit shall be 'scrsurv' class.")}
      if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}
      t1.formula<- fit$t1.formula
      t2.formula<- fit$t2.formula
      copulaparam<- fit$copulaparam
      if(indx[1]!=0){if(is.null(data))data<- fit$data}
    }
    if(indx[1]!=0){
      data<- prepData0(t1.formula,t2.formula, data)
      na.action<- data$na.action
      data<- data$data
      t1data <-  prepData(formula=t1.formula, data, na.action,bandwidth, discrete.level,xlimits)
      t2data <-  prepData(formula=t2.formula, data, na.action,bandwidth, discrete.level,xlimits)


    }else{
      t1data<- fit$t1data
      t2data<- fit$t2data

    }



    if(isTRUE(t1equal))stop("t1equal shall be FALSE for method 'np'.")

    if(length(all.vars(t1.formula[[3]]))>1)
    {stop("nonparametric estimation for stratified joint survival function is under construction.")}

    time1<- t1data$time
    event1<- t1data$event
    time2<- t2data$time
    event2<- t2data$event

    js0.t<- lapply(1:length(t1),function(k1)
      sapply(1:length(t2),function(k2)
        mean((time1>t1[k1])*(time2>t2[k2]))
        ))
    js0.s<- lapply(1:length(t1),function(k1)
      sapply(1:length(s2),function(k2)
        mean((time1>t1[k1])*(time2>s2[k2]))
      ))
    if(sum(1-event2)/n>0.05){
      Gfit<- survival::survfit(survival::Surv(time2,1-event2)~1)
      Gest.t2<- approx(Gfit$time,Gfit$surv,xout =  t2,method="constant",f = 0)$y
      Gest.s2<- approx(Gfit$time,Gfit$surv,xout =  s2,method="constant",f = 0)$y
      rm(Gfit)

      js.t<- lapply(js0.t,function(xx)ifelse(xx==0& Gest.t2==0,0, xx/ Gest.t2))
      js.s<- lapply(js0.s,function(xx)ifelse(xx==0& Gest.s2==0,0, xx/ Gest.s2))
    }else{
      js.t<- js0.t
      js.s<- js0.s
    }
    names(js.t)<- names(js.s)<- paste0("t1=",t1)
    condsurv<- lapply(1:length(t1), function(k)tcrossprod(js.t[[k]],1/js.s[[k]]))
    condsurv<- lapply(condsurv,t2s2matNA,t2=t2,s2=s2)

    condsurv<- lapply(1:length(t1), function(k)t1s2matNA(x=condsurv[[k]],tt1=t1[k],s2))
    names(condsurv)<- paste0("t1=",t1)
    return(list(csurv = condsurv,t1=t1,t2=t2,s2=s2))
  }
  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if(is.null(fit)){stop("fit shall be 'scrsurv' class.")}
    if(!inherits(fit,"scrsurv")){stop("fit shall be 'scrsurv' class.")}

    surv1<- fit$t1.surv
    surv2<- fit$t2.surv
    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam

    if (fit$t1.surv$covariate.type==1){
      s1.t1<- approx(x=c(0,surv1$time),y=c(1,surv1$surv),xout = t1,method = "constant",
                     yleft = 1,yright = min(surv1$surv),f = 0,na.rm = TRUE,ties = mean)$y

      s2.t2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = t2,method = "constant",
                     yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE,ties = mean)$y

      s2.s2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = s2,method = "constant",
                     yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE,ties = mean)$y



    }  else{

      strata.vars <- fit$t1.surv$discrete.predictors
      X.formula <- update(formula(fit$t1.formula),NULL~.)
      if (!all(match(all.vars(X.formula),names(newdata),nomatch=FALSE)))
        stop("Arg newdata does not contain all the covariates used for fitting. \n\nfitted variables: ", paste(all.vars(X.formula),collapse=", "),"\nnewdata contains:",ifelse(length(names(newdata))==0," nothing",names(newdata)))
      requested.X <- newdata[,all.vars(X.formula),drop=FALSE]


      fit.strata <- interaction(fit$stratas[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      requested.strata <- interaction(requested.X[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      fit.levels <- as.character(unique(fit.strata))
      NS <- length(fit.levels)
      fit.strata <- factor(fit.strata,levels=levels(fit.strata),labels=1:NS)
      requested.strata <- factor(requested.strata,levels=fit.levels,labels=1:NS)
      strata<- rep(fit.strata,times=fit$t1.surv$size.strata)
      copulaparam<- copulaparam[requested.strata]


      survmat1<- data.frame(time=surv1$time,surv=surv1$surv,strata=strata)
      survfn1<- lapply(fit.strata,function(k){
        sv<- survmat1[survmat1[,3]==k,]
        stats::stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
      })
      psurv1<- survfn1[requested.strata]
      s1.t1<- lapply(psurv1,function(x)x(t1))


      survmat2<- data.frame(time=surv2$time,surv=surv2$surv,strata=strata)
      survfn2<- lapply(fit.strata,function(k){
        sv<- survmat2[survmat2[,3]==k,]
        stats::stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
      })
      psurv2<- survfn2[requested.strata]
      s2.t2<- lapply(psurv2,function(x)x(t2))

      s2.s2<- lapply(psurv2,function(x)x(s2))

    }

    if(!isTRUE(t1equal)&fit$t1.surv$covariate.type==1){
      js.t<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.t2))


      js.s<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.s2))
    }
    if(isTRUE(t1equal)&fit$t1.surv$covariate.type==1){
      # js.t<-  lapply(1:length(t1),function(k)
      #   archmCopulaLink_dev(copulafam=copulafam,param=copulaparam, js=js.t[[k]]))
      js.t<-  lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param = copulaparam,p1 =s1.t1[k1],p2 =s2.t2,mode = "1"))

      # js.s<-  lapply(1:length(t1),function(k)
      #   archmCopulaLink_dev(copulafam=copulafam,param=copulaparam, js=js.s[[k]]))

      js.s<- lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param = copulaparam,p1 =s1.t1[k1],p2 =s2.s2,mode = "1"))

    }

    if(!isTRUE(t1equal)&fit$t1.surv$covariate.type!=1){
      js.t<- lapply(1:length(s1.t1),function(k)
      {lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param =copulaparam[k],p1 = s1.t1[[k]][k1],p2 = s2.t2[[k]]))
      })

      js.s<- lapply(1:length(s1.t1),function(k)
      {lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param =copulaparam[k],p1 = s1.t1[[k]][k1],p2 =  s2.s2[[k]]))
      })


    }
    if(isTRUE(t1equal)&fit$t1.surv$covariate.type!=1){
      js.t<- lapply(1:length(s1.t1),function(k)
      {lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param =copulaparam[k],p1 = s1.t1[[k]][k1],p2 = s2.t2[[k]],mode = "1"))
      })

      js.s<- lapply(1:length(s1.t1),function(k)
      {lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param =copulaparam[k],p1 = s1.t1[[k]][k1],
                   p2 =  s2.s2[[k]],mode = "1"))
      })


    }


    if(fit$t1.surv$covariate.type==1){
      names(js.t)<- names(js.s)<- paste0("t1=",t1)

      condsurv<- lapply(1:length(t1), function(k)tcrossprod(js.t[[k]],1/js.s[[k]]))
      condsurv<- lapply(condsurv,t2s2matNA,t2=t2,s2=s2)

      condsurv<- lapply(1:length(t1), function(k)t1s2matNA(x=condsurv[[k]],tt1=t1[k],s2))
      names(condsurv)<- paste0("t1=",t1)

    }else{

      condsurv<- lapply(1:length(js.t),function(x)
                        lapply(1:length(t1), function(k)tcrossprod(js.t[[x]][[k]],1/js.s[[x]][[k]])))
      condsurv<- lapply(1:length(js.t),function(x)lapply(condsurv[[x]],t2s2matNA,t2=t2,s2=s2))

      condsurv<- lapply(1:length(js.t),function(x){
        a<- lapply(1:length(t1), function(k)t1s2matNA(x=condsurv[[x]][[k]],tt1=t1[k],s2))
        names(a)<- paste0("t1=",t1)
        a})
      names(condsurv)<- paste(strata.vars,requested.X[,1],sep="=")
    }
    return(list(csurv = condsurv,t1=t1,t2=t2,s2=s2))
  }

}

#' @noRd
t2s2matNA<- function(x,t2,s2){

  colnames(x)<- paste0("s2=",s2)
  rownames(x)<- paste0("t2=",t2)


  tsind<- sapply(s2,function(xx)xx>t2)
  x[tsind]<- NA
  return(x)
}

#' @noRd
t1s2matNA<- function(x,tt1,s2){
  x[,tt1>s2]<- NA
  return(x)
}


