#################################################################
##   CopulaSCR R package by Tonghui Yu Copyright (C) 2024
##
## The below functions are used for association analysis for semi-competing
## risks data with single intermediate event time
#################################################################
#' @title Concordance estimator for association parameter
#' @description  Fitting a semiparametric copula-based model with pre-specified Archimedean copula.
#' The copula parameter is solved from a generalized concordance estimating equations proposed by Lakhal et al. (2008).
#' @aliases scrassonp
#' @usage scrassonp(t1.formula, t2.formula, data = parent.frame(),
#'      equalweight = FALSE, a = 0, b = 0,positive = TRUE, tol = 1e-5,
#'      copulafam = c("clayton", "frank", "joe", "gumbel", "amh"), se = FALSE,
#'      se.method = c("bootstrap", "resampling"), B = 0, seed = NULL)
#'
#' @param data an optional data.frame in which to interpret the variables occurring in the formula.
#' @param t1.formula a formula expression for the nonterminal event time,
#'     of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object with right censoring.
#'     See the documentation of \code{coxph} and \code{formula} for details.
#' @param t2.formula a formula expression for the terminal event time,
#'     of the form \code{response ~ predictors}.
## #' @param weights an optional vector of observation weights.
## #' @param na.action the na.action attribute, if any, that was returned by the na.action routine.
#' @param equalweight if \code{equalweight=TRUE} the weight function (see Details below) is \eqn{w = 1}.
#' Otherwise, numeric values in the arguments \code{a} and \code{b} will be used to compute the weight function.
#' @param copulafam a character string specifying the family of an Archimedean copula.
#' Currently supported families are "frank", "clayton", "amh", "gumbel", and "joe".
#' @param se Whether or not compute standard error (SE) for the concordance estimate
#' of copula parameter.
#' @param se.method the method used for computing SE for association parameter.
#' @param B a numeric value specifies the number of resampling/Bootstrap replicates.
#'     When B = 0, only the point estimate of copula parameter will be displayed.
## #' @param model If model= TRUE, the output includes \code{data}.
#' @param a positive constant to dampen $w(x, y) $ for large $x$ and $y$ (see Details below).
#' If \code{a = NULL} choose \code{quantile(T1,0.95)} for \eqn{a}
#' Default is \code{a=0}.
#' @param b positive constant dampen $w(x, y) $ for large $x$ and $y$ (see Details below).
#' If \code{b = NULL} choose \code{quantile(T2,0.95)} for \eqn{b}.
#' Default is \code{b=0}. When both \code{a=0} and \code{b=0}, the the weight function \eqn{w = 1},
#' namely, \code{equalweight=TRUE}.
#' @param tol the desired accuracy.
#' @param seed Integer. Specify seeds for the random generator to ensure
#' reproducibility of resampling/bootstrap replicates.
#' @param positive whether or not put positive constraint on the Kendall' tau parameter.
#'
#' @return An object of class "\code{scrassonp}" representing the fit.
#' The \code{scrassonp} object is a list containing at least the following components:\cr
#' \describe{
#' \item{Call}{An object of class \code{call}.}
#' \item{tau}{Estimated Kendall's tau corresponding to estimated copula parameter.}
#' \item{copulaparam}{Estimated copula parameter in the Archimedean copula.}
#' \item{copulafam}{Prespecified Archimedean copula.}
#' \item{zetas}{Resampling replicates.}
#' \item{param.boot}{Estimated copula parameter under each set of
#' resampling/Bootstrap replicates.}
#' \item{param.se}{SE for estimated copula parameter.}
#' \item{tau.boot}{Estimated Kendall's tau under each set of resampling/Bootstrap replicates.}
#' \item{tau.se}{SE for estimated Kendall's tau.}
#' \item{...}{The object will also contain the input arguments.}
#' }
#'
#' @details See details in Lakhal et al. (2008) and Yu et al. (2024)
#'
#' @references Fine, J. P., Jiang, H., and Chappell, R. (2001).
#' On semi-competing risks data. \emph{Biometrika}, \bold{88}(4):907–919.
#' @references Lakhal, L., Rivest, L.-P., and Abdous, B. (2008).
#' Estimating survival and association in a semicompeting risks model.
#' \emph{Biometrics}, \bold{64}(1):180–188.
#' @references Yu, Tonghui, Xiang, Liming, Chen, Chixiang, and Chiou, Sy Han (2024). CopulaSCR:  .
#' \emph{Working paper.}
## #' \emph{Journal of XXX}, \bold{105}(5): 1--34.
#'
#' @export
#' @examples
#' set.seed(12345)
#' simdata<- simSCR(n = 100,tau = 0.5, copulafam = "frank")
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
#'                     t2.formula=Surv(T2, event2) ~ 1,
#'                     data = simdata, copulafam="frank",
#'                     a=quantile(simdata$T1,0.9), b=quantile(simdata$T2,0.9),
#'                     B=20,seed = 12345,se = TRUE,se.method = "resampling")
#' cat("tau.est =", fitasso$tau, "tau.se = ",fitasso$tau.se)
#' # tau.est = 0.5216944 tau.se =  0.0891336
#'
#' set.seed(12345)
#' simdata<- simSCRtr(n = 100,K=3, tau = c(0.3,0.5,0.6), copulafam = "frank",
#'                    params=list(marginsDist = rep("exp",2),
#'                                rate1=c(0.5,1,1.2),rate2=1))
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ tr,
#'                     t2.formula=Surv(T2, event2) ~ tr,
#'                     data=simdata, copulafam="frank",
#'                     a=quantile(simdata$T1,0.9), b=quantile(simdata$T2,0.9),
#'                     B=20,seed = 12345,se = TRUE,se.method = "resampling")
#' cat("tau.est =", round(fitasso$tau,2), "tau.se = ",round(fitasso$tau.se,2))
#' # tau.est = 0.39 0.51 0.69 tau.se =  0.38 0.14 0.1
#' @seealso  \code{\link{scrsurv}}
#'
scrassonp<- function(t1.formula, t2.formula,data= parent.frame(),
                     equalweight = FALSE,a=0,b=0,positive=TRUE,tol=1e-5,
                     copulafam=c("clayton","frank","joe","gumbel","amh"),
                     se=FALSE,se.method = c("bootstrap","resampling"),
                     B=0,seed=NULL){

  Call <- match.call()
  xlimits <- 5
  discrete.level <- 20
  bandwidth <- NULL
  model<- TRUE
  if(B==0) {se= FALSE}
  if(B>0){se=TRUE}
  if(isTRUE(se)){
    se.method<- match.arg(se.method)
  }else{
    se.method<- NULL
  }

  copulafam<- match.arg(copulafam)

  indx <- match(c('t1.formula', 't2.formula','data',  'weights'),
                names(Call), nomatch=0)
  # if (indx[1]==0) stop("a data argument is required")
  if (indx[1]==0) stop("a t1.formula argument is required")
  if (indx[2]==0) stop("a t2.formula argument is required")

  newform<- prepformula(t1.formula,t2.formula)
  t1.formula<- newform$form1
  t2.formula<- newform$form2
  data<- prepData0(t1.formula,t2.formula, data)
  na.action<- data$na.action
  data<- data$data
  t1data <-  prepData(formula=t1.formula, data, na.action,bandwidth, discrete.level,xlimits)
  t2data <-  prepData(formula=t2.formula, data, na.action,bandwidth, discrete.level,xlimits)
  if(!identical(names(t2data$time),names(t1data$time))){
    stop("Please correct the formulas.")
  }else if(!identical(t2data$size.strata,t1data$size.strata)){
    stop("Please correct the formulas.")
  }else {
    out<- scrassonp0(t1data=t1data,t2data = t2data, equalweight = equalweight,
                     copulafam=copulafam,se = se,se.method = se.method,B=B,a = a,
                     b = b,tol = tol,seed = seed,positive = positive)
  }
  out$t1data<- t1data
  out$t2data<- t2data
  if(isTRUE(model)){
    out$data<- data

  }
  out$Call<- Call
  out$t1.formula<- t1.formula
  out$t2.formula<- t2.formula
  class(out)<- "scrassonp"
  return(out)

}

#' @aliases scrassonp0
#' @noRd
scrassonp0<- function(t1data,t2data, equalweight = FALSE,copulafam,
                      se=FALSE,se.method,B=100,a=0,b=0,tol=1e-5,
                      seed=NULL,positive=TRUE){

  time1<- t1data$time
  event1<- t1data$event
  time2<- t2data$time
  event2<- t2data$event

  n<- length(time1)

  casewt<- t1data$weights
  if (is.null(casewt)) {
    casewt <- rep(1.0, n)
  }  else {
    if (!is.numeric(casewt)) stop("weights must be numeric")
    if (any(!is.finite(casewt))) stop("weights must be finite")
    if (any(casewt <0)) stop("weights must be non-negative")
    casewt <- as.numeric(casewt)  # transform integer to numeric
  }


  size.strata<- t1data$size.strata
  nstrata<- t1data$nstrata
  if(nstrata==1){
    subsets<- NULL
  }else{
    sumstrata<- c(0,cumsum(size.strata))
    subsets<- lapply(1:nstrata,function(k) (sumstrata[k]+1):(sumstrata[k+1]))
  }


  if(nstrata>1){
    copulaparam<- c()
    for (k in 1:nstrata) {
      copulaparam<- c(copulaparam,
                      assonpfit( time1 = time1,event1 = event1,
                                 time2 = time2,event2 = event2,
                                 copulafam = copulafam,equalweight=equalweight,
                                 a=a,b= b,tol=tol,zeta=NULL,
                                 positive=positive,subset =  subsets[[k]]))

    }
    tau.est<- sapply(copulaparam,Caltau,copulafam = copulafam)
    names(tau.est)<- paste0("strata",1:nstrata)
  }else{
    copulaparam<- assonpfit( time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                             copulafam = copulafam,equalweight=equalweight,
                             a=a,b= b,tol=tol,zeta=NULL,positive=positive)
    tau.est<- Caltau(copulafam = copulafam,copulaparam = copulaparam)

  }



  out<- list(tau=tau.est,copulaparam=copulaparam,copulafam = copulafam)


  out$se<- se
  if(isTRUE(se)&!is.null(seed)&is.numeric(seed)){set.seed(seed)}

  if(isTRUE(se)){

    out<- scrassonp_se(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                       out = out,se.method = se.method,B=B,
                       equalweight =equalweight,a=a,b=b,tol=tol,
                       positive = positive,nstrata = nstrata,
                       subsets = subsets)
  }

  return(out)
}


#' @aliases assonpfit
#' @noRd
assonpfit<- function( time1,event1,time2,event2,equalweight,copulafam,
                      a,b,tol=1e-5,zeta=NULL,positive=TRUE,subset=NULL){

  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(zeta)){
    zeta<- zeta[subset]
  }
  n<- length(event2)
  if(is.null(zeta)){
    zeta2<- 1
    }else{
      zeta2<- utils::combn( zeta,m=2)
      zeta2<- zeta2[1,]*zeta2[2,]
    }

  t1.pair<- utils::combn(time1,m=2)
  t2.pair<- utils::combn(time2,m=2)
  event1.pair<- utils::combn(event1,m=2)
  event2.pair<- utils::combn(event2,m=2)

  t1.pairmin<- apply(t1.pair,2,min)
  t2.pairmin<- apply(t2.pair,2,min)

  pair<- ((t1.pair[1,]>t1.pair[2,])*event1.pair[2,]+
            (t1.pair[1,]<t1.pair[2,])*event1.pair[1,])
  pair<- pair*((t2.pair[1,]<t2.pair[2,]) *event2.pair[1,]
               +(t2.pair[1,]>t2.pair[2,])*event2.pair[2,])


  Delta<- (t1.pair[1,]-t1.pair[2,])*(t2.pair[1,]-t2.pair[2,])
  Delta<- 1* (Delta>0 )


  if(equalweight){
    w<- 1
  }else{
    if(is.null(a)){a<- stats::quantile(time1,0.95)}
    if(is.null(b)){b<- stats::quantile(time2,0.95)}
    w_inv<- sapply(1:length(t1.pairmin),function(kk)
      mean((time1>=min(a,t1.pairmin[kk]))*(time2>=min(b,t2.pairmin[kk]))))
    w<- 1/w_inv
    rm(w_inv)
  }

  w<- w*zeta2

  if(copulafam=="clayton"){
    numer<- sum(w*Delta *pair,na.rm = TRUE)
    denom<- sum(w*(1-Delta) *pair,na.rm = TRUE)
    copulaparam<- numer/denom -1

  }else{
    js0<- sapply(1:length(t1.pairmin),function(kk)
      mean((time1>t1.pairmin[kk])*(time2>t2.pairmin[kk])))
    if(sum(1-event2)/n>0.05){
      Gfit<- survival::survfit(survival::Surv(time2,1-event2)~1,weights = zeta)
      Gest<- stats::approx(c(0,Gfit$time),c(1,Gfit$surv),yright = min(Gfit$surv),ties = mean,
                    xout =  t2.pairmin,method="constant",f = 0)$y
      # Gfit<- prodlim::prodlim(Surv(time2,event2)~1,reverse = TRUE,caseweights = zeta)
      # Gest<-  predict(Gfit,times = t2.pairmin)
      rm(Gfit)

      js<- js0/ Gest
      js[js0==0& Gest==0]<- 0
    }else{
      js<- js0
    }

    if(isTRUE(positive)){
      paramin<-copulabd_pos(copulafam)
    }else{
      paramin<- copulabd(copulafam)
    }


    # copulaparam<- stats::uniroot(f = function(par){
    #   scrnpGEE(gamma=par,copulafam = copulafam,w = w,js = js,pair = pair,
    #            Delta = Delta)},
    #   interval = paramin,trace = FALSE,maxiter =maxiter,tol = tol)$root

    copulaparam<- stats::optimize(f = function(par){
      abs(scrnpGEE(gamma=par,copulafam = copulafam,w = w,js = js,pair = pair,
                   Delta = Delta))},
      interval = paramin,tol = tol)$minimum


  }
  return(copulaparam)

}


#' @aliases scrassonp_se
#' @noRd
scrassonp_se<- function(time1,event1,time2,event2,out,se.method,B=100,
                        equalweight,a,b,tol,positive,nstrata,subsets){
  n<- length(time1)
  out$se.method<- se.method
  out$B<- B
  copulafam<- out$copulafam
  copulaparam.boot<- c()
  cat("Please wait for a while...")
  if(nstrata>1){

    if(se.method =="resampling"){
      zetas<- c()
      for(j in 1:B){
        zeta<-rexp(n,rate=1)
        params<- c()
        for (k in 1:nstrata) {
          params<- c(params,
                     assonpfit( time1 = time1,event1 = event1,
                                time2 = time2,event2 = event2,
                                copulafam = copulafam,equalweight=equalweight,
                                a=a,b= b,tol=tol,zeta=zeta,positive=positive,
                                subset =  subsets[[k]]))


        }
        copulaparam.boot<- rbind(copulaparam.boot,params)
        zetas<- rbind(zetas,zeta)
      }
      out$zetas<- zetas
    }else{
      bootind<- c()
      for(j in 1:B){

        params<- c()
        name.indd<- c()
        for (k in 1:nstrata) {
          subset <-  subsets[[k]]
          indd<- sort(sample(length(subset),replace = TRUE))

          params<- c(params,
                     assonpfit(time1 = time1[subset][indd],event1 = event1[subset][indd],
                               time2 = time2[subset][indd],event2 = event2[subset][indd],
                               copulafam = copulafam,equalweight=equalweight,
                               a=a,b= b,tol=tol,zeta=NULL,positive=positive))

          name.indd<- c(name.indd,names(time1)[subset][indd])
        }
        copulaparam.boot<- rbind(copulaparam.boot,params)
        bootind<- rbind(bootind,name.indd)
      }
      out$bootind<- bootind
    }
    colnames(copulaparam.boot)<- paste0("strata",1:nstrata)

    out$param.boot<- copulaparam.boot
    out$param.se<- apply(copulaparam.boot,2,mad)
    out$tau.boot<- apply(copulaparam.boot,2,
                         function(x)sapply(x,Caltau,copulafam =copulafam))
    out$tau.se<- apply(out$tau.boot,2,mad)

  }else{
    if(se.method =="resampling"){
      zetas<- c()
      for(k in 1:B){
        zeta<-rexp(n,rate=1)
        copulaparam.boot<- c(copulaparam.boot,
                             assonpfit( time1 = time1,event1 = event1,
                                        time2 = time2,event2 = event2,
                                        copulafam = copulafam,equalweight=equalweight,
                                        a=a,b= b,tol=tol,zeta=zeta,positive=positive))
        zetas<- rbind(zetas,zeta)

      }
      out$zetas<- zetas

    }else{
      bootind<- c()
      for(k in 1:B){
        indd<- sort(sample(1:n,replace = TRUE))

        copulaparam.boot<- c(copulaparam.boot,
                             assonpfit( time1 = time1[indd],event1 = event1[indd],
                                        time2 = time2[indd],event2 = event2[indd],
                                        copulafam = copulafam,equalweight=equalweight,
                                        a=a,b= b,tol=tol,zeta=NULL,positive=positive))

        bootind<- rbind(bootind,indd)
      }
      out$bootind<- bootind
    }

    out$param.boot<- copulaparam.boot
    out$param.se<- sd(copulaparam.boot)
    out$tau.boot<- sapply(copulaparam.boot, Caltau,copulafam =copulafam)
    out$tau.se<- sd(out$tau.boot)
  }




  return(out)
}


#' @aliases scrnpGEE
#' @noRd
scrnpGEE<- function(gamma,copulafam,w,js,pair,Delta){

  chr<- npCHRfit(s=js,param = gamma,copulafam = copulafam)
  mean(w *pair* (Delta-  chr/(chr+1)),na.rm = TRUE)

  # mean((w *(Delta-  chr/(chr+1)))[pair==1],na.rm = TRUE)
}

#' @aliases npCHRfit
#' @noRd
npCHRfit<- function(s,param,copulafam){
  if(copulafam=="amh"){
    phi<- (2*param*s +1-param)/(param*s +1-param)
  }
  if(copulafam=="clayton"){
    phi<- param +1
  }
  if(copulafam=="frank"){
    phi<- (param*s)/(1-exp(-param*s))
    phi[s==0]<- 1
  }
  if(copulafam=="gumbel"){
    phi<- 1 + (1-param)/log(s)
  }
  if(copulafam=="joe"){
    phi<-  s*(param-1+(1-s)^param)/(1-s) /(1-(1-s)^param)
    phi[s==0]<- 1
  }
  return(phi)
}

#' @aliases copulabd_pos
#' @noRd
copulabd_pos<- function(copulafam){
  lb<- switch(copulafam,
              amh = -1,
              clayton = 0,
              frank = 0.0001,
              gumbel = 1,
              joe = 1)
  ub<- switch(copulafam,
              amh = 1,
              clayton = 1000,
              frank = 1000,
              gumbel = 1000,
              joe = 1000)
  return(c(lb,ub))
}


