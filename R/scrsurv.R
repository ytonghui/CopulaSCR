#' @title Nonparametric estimation of marginal survival functions
#'
#' @description Fitting a semiparametric copula-based model with pre-specified
#'   Archimedean copula. Plugging estimated copula parameter solving from a
#'   generalized concordance estimating equations proposed by Lakhal et al.
#'   (2008), the marginal survival function of nonterminal event time can be
#'   estimated using one of three available methods: "FJC" for Fine et al.
#'   (2001)'s estimator, "JFKC" for Jiang et al. (2005)'s estimator, and "LRA"
#'   for Lakhal et al. (2008)'s estimator. Or alternatively, a Kernal-based
#'   estimator ("Ker") is also provided. The marginal survival function of
#'   terminal event time can be estimated using one of two available methods in:
#'   "JFKC" for Jiang et al. (2005)'s estimator, and the classical Kaplan-Meier
#'   estimator (see survival package).
#'
#' @aliases scrsurv
#' @usage scrsurv(fit, t1.formula, t2.formula, data = parent.frame(), tau,
#'   copulafam = c("clayton", "frank", "joe", "gumbel", "amh"),
#'   method = c("JFKC", "LRA", "FJC", "Ker"), surv2km = TRUE,
#'   B = 100, se.method = c("bootstrap","resampling"), conf.int = TRUE,
#'   conftype = 1, seed = NULL)
#'
#' @param data an optional data frame in which to interpret the variables
#'   occurring in the formula.
#' @param fit an optional object of class  "scrassonp".
#' @param t1.formula  an optional formula expression for the nonterminal event
#'   time, of the form \code{response ~ predictors}. The \code{response} is a
#'   \code{Surv} object with right censoring. See the documentation of
#'   \code{coxph} and \code{formula} for details.
#' @param t2.formula  an optional formula expression for the terminal event
#'   time, of the form \code{response ~ predictors}.
## #' @param weights an optional vector of observation weights.
## #' @param na.action the na.action attribute, if any, that was returned by the
#'   na.action routine.
#' @param tau an optional numeric value of Kendall's tau parameter in the
#'   Archimedean copula.
#' @param copulafam a character string specifying the family of an Archimedean
#'   copula. Currently supported families are "frank", "clayton", "amh",
#'   "gumbel", and "joe".
#' @param method a character string specifying the method for computing survival
#'   curves. Currently supported estimators for nonterminal event time's
#'   survival function are "JFKC", "LRA", "FJC", "Ker", available estimators for
#'   terminal event time's survival function are "JFKC" estimator (if method =
#'   "JFKC" and surv2km = FALSE), and the Kaplan-Meier estimator (if method =
#'   "LRA", "FJC", "Ker" or surv2km = TRUE)
#' @param surv2km whether or not applying the Kaplan-Meier method for estimating
#' survival curve of terminal event time. surv2km = TRUE when method != "JFKC".
#' @param se.method The method used for computing SE for marginal survival curves.
#' @param B  A numeric value specifies the number of resampling/Bootstrap replicates.
#'     When B = 0, only the point estimate of copula parameter will be displayed.
#' @param conf.int Whether or not compute 95\% confidence interval for marginals.
#' @param conftype The type of confidence interval can be specified with values
#'   1, 2, 3, and 4.
#'
#'   \code{type = 1} and \code{type = 3} indicate Wald-type confidence intervals
#'   derived from resampling/Bootstrap method. \code{type = 2} and \code{type =
#'   4} correspond to  quantile-based confidence intervals from the
#'   resampling/bootstrap method. \code{type = 1} and \code{type = 2} indicate
#'   that the copula parameter is prespecified. \code{type = 3} and \code{type =
#'   4} indicate the use of resampling/Bootstrap method in the estimation of
#'   both the copula parameter and the marginals.
## #' @param model If model= TRUE, the output includes the input data.
#' @param seed Integer. Specify seeds for the random generator to ensure
#' reproducibility of resampling/bootstrap replicates.
#'
#'
#'
#' @return An object of class "\code{scrsurv}" representing the fit.
#' The \code{scrsurv} object is a list containing at least the following components:\cr
#' \describe{
#' \item{t1.surv}{A list containing "time", "n.risk", "n.event", "n.censor",
#' "surv", "cumhaz", "type", "std.err", "std.chaz", "lower", "upper"
#' (See \code{survfit.object()} in \code{survival} package).}
#' \item{t2.surv}{An object of class "list" if surv2km = FALSE or class "survfit"
#' if surv2km = TRUE (See \code{survfit.object()} in \code{survival} package).}
#' \item{copulaparam}{Copula parameter in the Archimedean copula.}
#' \item{copulafam}{Prespecified Archimedean copula structure.}
#' \item{tau}{Kendall's tau corresponding to copula parameter.}
#' \item{Call}{An object of class \code{call}.}
#' \item{copulaparams}{Estimated copula parameter under each set of resampling/Bootstrap replicates.}
#' \item{t1data}{The data for nonterminal event times of class "list".
#' It includes observed times (\code{"time"}), censoring indicators (\code{"event"}),
#' covariates (\code{"X"})}
#' \item{t2data}{The data for terminal event times of class "list". }
#' \item{...}{}
#' }
#'
#' @export
#'
#' @references Fine, J. P., Jiang, H., and Chappell, R. (2001). On
#'   semi-competing risks data. \emph{Biometrika}, \bold{88}(4):907–919.
#' @references Jiang, H., Fine, J. P., Kosorok, M. R., and Chappell, R. (2005).
#'   Pseudo self-consistent estimation of a copula model with informative
#'   censoring. \emph{Scandinavian Journal of Statistics}, \bold{32}(1):1–20.
#' @references Lakhal, L., Rivest, L.-P., and Abdous, B. (2008). Estimating
#'   survival and association in a semicompeting risks model. \emph{Biometrics},
#'   \bold{64}(1):180–188.
#' @references Nevo D, Gorfine M (2022). Causal inference for semi-competing
#'   risks data. \emph{Biostatistics}, \bold{23}(4), 1115–1132.
#' @references Yu, Tonghui, Xiang, Liming, Chen, Chixiang, and Chiou, Sy Han
#'   (2024). CopulaSCR:  .
#'   \emph{Working paper.}
###'   \emph{Journal of XXX}, \bold{105}(5): 1--34.
#'
#'
#'
#' @seealso  \code{\link{scrassonp}}
#' @importFrom survival survfit Surv
#'
#' @examples
#' set.seed(12345)
#' simdata<- simSCR(n = 50,tau = 0.5, copulafam = "frank")
#' fitasso<- scrassonp(simdata, t1.formula=Surv(T1, event1) ~ 1,
#'          t2.formula=Surv(T2, event2) ~ 1,copulafam="frank",
#'          a=quantile(simdata$T1,0.9),
#'          b=quantile(simdata$T2,0.9),B=10,seed = 12345,
#'          se = TRUE,se.method = "resampling")
#'
#' \dontrun{
#' cat("tau.est =", fitasso$tau, "tau.se = ",fitasso$tau.se)
#' fit21<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=10,
#' se.method = "resampling",conf.int=TRUE,conftype =1,seed = 12345)
#' plot(fit21,conf.int = TRUE)
#' fit22<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=fitasso$B,
#' se.method = fitasso$se.method,conf.int=TRUE,conftype =3)
#' plot(fit22,conf.int = TRUE)
#'
#' #####  estimation of nonparmatric marginals with specified Kendall's tau
#' fit1<- scrsurv(t1.formula=Surv(T1, event1) ~ 1, t2.formula=Surv(T2, event2) ~ 1,
#' data =simdata, tau=0.5,copulafam="frank", B=10,
#' method= "FJC",se.method = "bootstrap",conf.int=TRUE,seed = 12345)
#' plot(fit1,conf.int = TRUE)
#'
#'
#' set.seed(12345)
#' simdata2<- simSCRtr(n = 500,K=3, tau = c(0.3,0.5,0.6), copulafam = "frank",
#'                     params=list(marginsDist = rep("exp",2), rate1=c(0.5,1,1.2),rate2=1))
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ tr,
#'                     t2.formula=Surv(T2, event2) ~ tr,
#'                     data=simdata2, copulafam="frank",
#'                     a=quantile(simdata2$T1,0.9), b=quantile(simdata2$T2,0.9),
#'                     B=10,seed = 12345,se = TRUE,se.method = "resampling")
#' cat("tau.est =", round(fitasso$tau,2), "tau.se = ",round(fitasso$tau.se,2))
#'
#' fit21tr<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=10,seed = 12345,
#'                   se.method = "resampling",conf.int=TRUE,conftype =1)
#' print(fit21tr)
#' plot(fit21tr,conf.int=TRUE)
#' }



scrsurv <- function(fit, t1.formula, t2.formula, data = parent.frame(),
                    tau,copulafam=c("clayton","frank","joe","gumbel","amh"),
                    method= c("JFKC","LRA","FJC","Ker"),surv2km=TRUE,B=100,
                    se.method = c("bootstrap","resampling"),
                    conf.int=TRUE, conftype=1,seed=NULL ) {
  Call <- match.call()
  bandwidth = NULL
  model = TRUE
  xlimits <- 5
  discrete.level <- 20
  method<- match.arg(method)
  method<- substr(toupper(method),1,1)
  copulafam<- match.arg(copulafam)

  if(method=="K"){tau<- NULL}
  if(method!="J"){surv2km<- TRUE}

  if(B==0) {conf.int<- FALSE}

  if (isTRUE(conf.int)) conf.int <- 0.95

  if (is.numeric(conf.int)) {
    if (!(0 < conf.int && conf.int < 1))
      conf.int <- 0.95
  }



  if(is.numeric(conf.int)){
    se.method<- match.arg(se.method)
    if(method!="J"){se.method <- "bootstrap"}
    if(se.method !="resampling"){se.method <- "bootstrap"}
  }else{
    se.method<- NULL
  }

  indx <- match(c('data','fit' ,'t1.formula', 't2.formula'), names(Call), nomatch=0)

  if (indx[2]==0) {
    if (indx[1]==0) stop("a data argument is required")
    if (indx[3]==0) stop("a t1.formula argument is required")
    if (indx[4]==0) stop("a t2.formula argument is required")
  }


  if (indx[2]!=0) {
    if(!inherits(fit,"scrassonp")){
      stop("fit shall be 'scrassonp' class.")
    }else if (is.null(fit$data)&indx[1]==0) {
      stop("a data argument is required")
    }else if (!is.null(fit$data)&indx[1]==0){
      data<- fit$data
    }

    if (!is.null(fit$t1.formula))  t1.formula<- fit$t1.formula
    if (!is.null(fit$t2.formula)) t2.formula<- fit$t2.formula

    copulafam<- fit$copulafam
    tau<- fit$tau
    copulaparam<- fit$copulaparam

    if(!is.null(fit$t1data)&!is.null(fit$t2data)){
      t1data <- fit$t1data
      t2data <- fit$t2data
    }else{
      newform<- prepformula(t1.formula,t2.formula)
      t1.formula<- newform$form1
      t2.formula<- newform$form2
      data<- prepData0(t1.formula,t2.formula, data)
      na.action<- data$na.action
      data<- data$data
      t1data <-  prepData(formula=t1.formula, data, na.action,bandwidth, discrete.level,xlimits)
      t2data <-  prepData(formula=t2.formula, data, na.action,bandwidth, discrete.level,xlimits)

    }

  }
  if (!(is.matrix(data) | is.data.frame(data))) {
    stop("data object must be a matrix or a data.frame")
  }

  if(indx[2]==0){
    newform<- prepformula(t1.formula,t2.formula)
    t1.formula<- newform$form1
    t2.formula<- newform$form2
    data<- prepData0(t1.formula,t2.formula, data)
    na.action<- data$na.action
    data<- data$data
    t1data <-  prepData(formula=t1.formula, data, na.action,bandwidth, discrete.level,xlimits)
    t2data <-  prepData(formula=t2.formula, data, na.action,bandwidth, discrete.level,xlimits)
  }


  if(is.numeric(conf.int)&!is.null(seed)&is.numeric(seed)){set.seed(seed)}

  if(!identical(names(t2data$time),names(t1data$time))){
    stop("Please correct the formulas.")
  }else if(!identical(t2data$size.strata,t1data$size.strata)){
    stop("Please correct the formulas.")
  }else {
    if(indx[2]==0){
      sout<- scrsurv0(t1data=t1data,t2data = t2data,method=method,
                      tau=tau,copulafam=copulafam,conf.int=conf.int,
                      se.method=se.method,B=B,conftype=conftype,
                      model=model,surv2km=surv2km)
    }else{
      sout<- scrsurv0(fit,t1data=t1data,t2data = t2data,method=method,
                      tau=tau,copulafam=copulafam,conf.int=conf.int,
                      se.method=se.method,B=B,conftype=conftype,
                      model=model,surv2km=surv2km)
    }


  }


  sout$t1.formula<- t1.formula
  sout$t2.formula<- t2.formula
  if(is.numeric(conf.int)){sout$se.method<- se.method  }
  if(is.numeric(conf.int)&conftype%in% 3:4){
    sout$se.method<- se.method
    sout$copulaparams = fit$param.boot
    sout$conftype<- conftype
    # if(se.method == "bootstrap") {
    #   sout$bootind = fit$bootind
    # }
    # if(se.method == "resampling") {
    #   sout$zetas = fit$zetas
    # }

  }

  if(isTRUE(model)){
    sout$t1data<- t1data
    sout$t2data<- t2data
    sout$data<- data
  }
  sout<-c( sout, list(na.action= t2data$na.action,Call=Call))
  class(sout)<- "scrsurv"
  return(sout)

}

#' @aliases scrsurv0
#' @noRd
scrsurv0<- function(fit,t1data,t2data,method,tau,copulafam,conf.int,se.method,B,
                    conftype,model,surv2km) {

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

  if(length(tau)==1&nstrata>1){
    tau<- rep(tau,nstrata)
  }
  if(method!="K"){
    copulaparam<- sapply(tau,Calitau,copulafam = copulafam)
  }else{
    copulaparam<- NULL
  }

  if(!is.numeric(conf.int)|conftype%in% 1:2){
    sout<- scrsurvfit(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                      t1data = t1data,t2data = t2data,
                      copulaparam = copulaparam,copulafam = copulafam,
                      model = model,casewt=casewt,surv2km=surv2km,
                      method=method,conf.int =conf.int, se.method=se.method,
                      B=B,conftype=conftype,size.strata=size.strata,nstrata=nstrata)
  }
  if(is.numeric(conf.int)&conftype%in% 3:4){
    if (!isTRUE(fit$se)){
      stop("a se argument in fit shall be TRUE")
    }
    if (fit$se.method!= se.method ){
      se.method<- fit$se.method
      # stop("se.method arguments in fit and scrsurv function shall be consistent.")
    }

    if (fit$B!= B ){
      # stop("'B' in 'CopulaSCRnp' and scrsurv function shall be consistent.")
      B<- fit$B
    }

    sout<- scrsurvfit(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                      t1data = t1data,t2data = t2data,
                      copulaparam = copulaparam,copulafam = copulafam,
                      copulaparams = fit$param.boot,bootind = fit$bootind,zetas = fit$zetas,
                      model = model,casewt=casewt,surv2km=surv2km,
                      method=method,conf.int =conf.int, se.method=se.method,
                      B=B,conftype=conftype,size.strata=size.strata,nstrata=nstrata)
  }

  if(nstrata>1){sout$stratas<- t1data$stratas}
  if(method!="K"){
    sout$copulaparam<- copulaparam
    sout$copulafam<- copulafam
    sout$tau<- tau
  }
  return(sout)
}

#' @aliases scrsurvfit
#' @noRd
scrsurvfit<- function(time1,event1,time2,event2,t1data,t2data,copulaparam,copulafam,
                      copulaparams = NULL,bootind =NULL,zetas =NULL,
                      model,casewt,surv2km,method,conf.int, se.method,B,
                      conftype,size.strata,nstrata){

  n<- length(time1)
  if(nstrata==1){
    subsets<- NULL
  }else{
    sumstrata<- c(0,cumsum(size.strata))
    subsets<- lapply(1:nstrata,function(k) (sumstrata[k]+1):(sumstrata[k+1]))
    rm(sumstrata)
    zetas<- NULL
    bootind<- NULL
  }

  if(is.numeric(conf.int)){
    zstar<- qnorm(((1 - conf.int))/2,lower.tail = FALSE)
  }
  if(method=="J"){
    if(nstrata ==1 ){
      sout<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                             copulaparam = copulaparam,copulafam = copulafam,
                             model = model,casewt=casewt,surv2km=surv2km)
      names(sout)<- c("t1.surv","t2.surv")
      if(is.numeric(conf.int)){
      cat("Please wait for a while...")
      sout<- scrsurvse_J(sout=sout,conftype = conftype,se.method = se.method,
                         B = B,time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                         copulaparam =copulaparam,copulaparams = copulaparams,
                         bootind = bootind,copulafam = copulafam,model = model,
                         casewt = casewt,zetas = zetas,
                         surv2km = surv2km,zstar = zstar)
      }
    }else{
      sout<- list()
      if(is.numeric(conf.int)){
        cat("Please wait for a while...")
      }
      for (k in 1:nstrata) {
        sout[[k]]<- scrsurvfit.JFKC(time1 = time1,event1 = event1,
                                    time2 = time2,event2 = event2,
                                    copulaparam = copulaparam[k],copulafam = copulafam,
                                    model = model,casewt=casewt,
                                    surv2km=surv2km,subset = subsets[[k]])
        names(sout[[k]])<- c("t1.surv","t2.surv")

        if(is.numeric(conf.int)){

          sout[[k]]<- scrsurvse_J(sout=sout[[k]],conftype = conftype,se.method = se.method,
                                  B = B,time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                  subset = subsets[[k]],
                                  copulaparam =copulaparam[k],copulaparams = copulaparams[,k],
                                  bootind = bootind[,k],copulafam = copulafam,model = model,
                                  casewt = casewt,zetas = zetas[[k]],
                                  surv2km = surv2km,zstar = zstar)
        }

      }

      sout<- surv_transform(survlist=sout,inverse=FALSE,conf.int = conf.int,
                            t1data = t1data,t2data = t2data)

    }


  }


  if(method!="J"){

    if(nstrata ==1 ){
      S2fitkm<-survival::survfit(formula = survival::Surv(time2,event2)~1,
                        model = model)


      S1fit<- switch (method,
                      "L" = scrsurvfit.LRA(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                           copulaparam = copulaparam,copulafam = copulafam,
                                           model = model,casewt=casewt),
                      "F" = scrsurvfit.FJC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                           copulaparam = copulaparam,copulafam = copulafam,
                                           model = model,casewt=casewt,S2=S2fitkm),
                      "K" = scrsurvfit.kernel(time1 = time1,event1 = event1,time2 = time2,
                                              event2 = event2,S2=S2fitkm))
      if(is.numeric(conf.int)){
        cat("Please wait for a while...")
        sout<- scrsurvse_noJ(S1fit = S1fit,S2fitkm = S2fitkm,method = method,
                             conftype = conftype,se.method = se.method,B = B,
                             time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                             copulaparam= copulaparam,copulaparams = copulaparams,
                             bootind = bootind,copulafam = copulafam,model = model,
                             casewt = casewt,zetas = zetas,zstar = zstar)

      }else{
        sout<- list(t1.surv=S1fit,t2.surv=S2fitkm)
      }


    }else{
      if(is.numeric(conf.int)){
        cat("Please wait for a while...")
      }
      sout<- S2fitkm<- S1fit<- list()
      for (k in 1:nstrata) {
        S2fitkm[[k]]<- survival::survfit(formula = survival::Surv(time2,event2)~1,
                                         model = model,subset = subsets[[k]])

        S1fit[[k]]<- switch (method,
                        "L" = scrsurvfit.LRA(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                             copulaparam = copulaparam[k],copulafam = copulafam,
                                             model = model,casewt=casewt,subset = subsets[[k]]),
                        "F" = scrsurvfit.FJC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                                             copulaparam = copulaparam[k],copulafam = copulafam,
                                             model = model,casewt=casewt,S2=S2fitkm[[k]],subset = subsets[[k]]),
                        "K" = scrsurvfit.kernel(time1 = time1,event1 = event1,time2 = time2,
                                                event2 = event2,S2=S2fitkm[[k]],subset = subsets[[k]]))

        if(is.numeric(conf.int)){
          sout[[k]]<- scrsurvse_noJ(S1fit = S1fit[[k]],S2fitkm = S2fitkm[[k]],method = method,
                               conftype = conftype,se.method = se.method,B = B,
                               time1 = time1,event1 = event1,time2 = time2,
                               event2 = event2,subset= subsets[[k]],
                               copulaparam= copulaparam[k],copulaparams = copulaparams[,k],
                               bootind = bootind[,k],copulafam = copulafam,model = model,
                               casewt = casewt,zetas = zetas[[k]],zstar = zstar)

        }else{
          sout[[k]]<- list(t1.surv=S1fit[[k]],t2.surv=S2fitkm[[k]])
        }
      }

      sout<- surv_transform(survlist=sout,inverse=FALSE,conf.int = conf.int,
                            t1data = t1data,t2data = t2data)

    }

  }

  sout$t1.surv$model.response<- t1data$model.response
  sout$t1.surv$model.matrix<- t1data[["model.matrix"]]
  sout$t1.surv$formula<- t1data$formula
  sout$t1.surv$covariate.type<- t1data$covariate.type


  sout$t2.surv$model.response<- t2data$model.response
  sout$t2.surv$model.matrix<- t2data[["model.matrix"]]
  sout$t2.surv$formula<- t2data$formula
  sout$t2.surv$covariate.type<- t2data$covariate.type
  return(sout)
}

#' @aliases scrsurvse_J
#' @noRd
scrsurvse_J<- function(sout,conftype,se.method,B,time1,event1,time2,event2,subset=NULL,
                       copulaparam,copulaparams=NULL,bootind=NULL,copulafam,model,
                       casewt,zetas=NULL,surv2km,zstar){
  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(casewt)){
    casewt<- casewt[subset]
  }

  S1list<- list()
  H1list<- list()
  S2list<- list()
  H2list<- list()

  n<- length(time1)
  if(conftype %in%1:2){
    if(se.method =="resampling"){

      for(k in 1:B){
        zeta<-rexp(n,rate=1)
        fit<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                              copulaparam = copulaparam,copulafam = copulafam,
                              model = model,casewt=casewt*zeta,surv2km=surv2km)

        S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
        H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
        if(!isTRUE(surv2km)){
          S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
          H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
        }
        if(k%%20==0){cat(round(k/B,1)*100,"%.")}
      }

    }else{
      for(k in 1:B){
        indd<- sort(sample(1:n,replace = TRUE))
        fit<- scrsurvfit.JFKC(time1 = time1[indd],event1 = event1[indd],
                              time2 = time2[indd],event2 = event2[indd],
                              copulaparam = copulaparam,copulafam = copulafam,
                              model = model,casewt=casewt[indd],surv2km=surv2km)

        S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
        H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
        if(!isTRUE(surv2km)){
          S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
          H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
        }
        if(k%%20==0){cat(round(k/B,1)*100,"%.")}
      }
    }

    S1list<- do.call(rbind,S1list)
    H1list<- do.call(rbind,H1list)
    sout$t1.surv$std.err<- apply(S1list,2,sd,na.rm=TRUE)
    sout$t1.surv$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
    if(conftype==2){
      sout$t1.surv$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
      sout$t1.surv$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
    }else if(conftype==1){
      sout$t1.surv$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                         zstar*sout$t1.surv$std.err,0))
      sout$t1.surv$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                         zstar*sout$t1.surv$std.err,0))
    }

    if(!isTRUE(surv2km)){
      S2list<- do.call(rbind,S2list)
      H2list<- do.call(rbind,H2list)
      sout$t2.surv$std.err<- apply(S2list,2,sd,na.rm=TRUE)
      sout$t2.surv$std.chaz<- apply(H2list,2,sd,na.rm=TRUE)
      if(conftype==2){
        sout$t2.surv$lower<- apply(S2list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        sout$t2.surv$upper<- apply(S2list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==1){
        sout$t2.surv$lower<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)-
                                           zstar*sout$t2.surv$std.err,0))
        sout$t2.surv$upper<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)+
                                           zstar*sout$t2.surv$std.err,0))
      }

    }
  }
  if(conftype %in%3:4){
    if(se.method =="resampling"){

      for(k in 1:B){
        if(is.null(zetas)){
          zeta<- rexp(n,rate=1)
        }else{
          zeta<- zetas[k,]
          }

        fit<- scrsurvfit.JFKC(time1 = time1,event1 = event1,time2 = time2,event2 = event2,
                              copulaparam = copulaparams[k],copulafam = copulafam,
                              model = model,casewt=casewt*zeta,surv2km=surv2km)

        S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
        H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
        if(!isTRUE(surv2km)){
          S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
          H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
        }
        if(k%%20==0){cat(round(k/B,1)*100,"%.")}
      }

    }else{
      for(k in 1:B){
        if(is.null(bootind)){
          indd<- sort(sample(1:n,replace = TRUE))
        }else{
          indd<- bootind[k,]
        }

        fit<- scrsurvfit.JFKC(time1 = time1[indd],event1 = event1[indd],
                              time2 = time2[indd],event2 = event2[indd],
                              copulaparam = copulaparams[k],copulafam = copulafam,
                              model = model,casewt=casewt[indd],surv2km=surv2km)

        S1list[[k]]<- stepfun(fit$s1$time,c(1,fit$s1$surv),right = FALSE)(sout$t1.surv$time)
        H1list[[k]]<- stepfun(fit$s1$time,c(0,fit$s1$cumhaz),right = FALSE)(sout$t1.surv$time)
        if(!isTRUE(surv2km)){
          S2list[[k]]<- stepfun(fit$s2$time,c(1,fit$s2$surv),right = FALSE)(sout$t2.surv$time)
          H2list[[k]]<- stepfun(fit$s2$time,c(0,fit$s2$cumhaz),right = FALSE)(sout$t2.surv$time)
        }
        if(k%%20==0){cat(round(k/B,1)*100,"%.")}
      }
    }

    S1list<- do.call(rbind,S1list)
    H1list<- do.call(rbind,H1list)
    sout$t1.surv$std.err<- apply(S1list,2,sd,na.rm=TRUE)
    sout$t1.surv$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
    if(conftype==4){
      sout$t1.surv$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
      sout$t1.surv$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
    }else if(conftype==3){
      sout$t1.surv$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                         zstar*sout$t1.surv$std.err,0))
      sout$t1.surv$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                         zstar*sout$t1.surv$std.err,0))
    }

    if(!isTRUE(surv2km)){
      S2list<- do.call(rbind,S2list)
      H2list<- do.call(rbind,H2list)
      sout$t2.surv$std.err<- apply(S2list,2,sd,na.rm=TRUE)
      sout$t2.surv$std.chaz<- apply(H2list,2,sd,na.rm=TRUE)
      if(conftype==4){
        sout$t2.surv$lower<- apply(S2list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        sout$t2.surv$upper<- apply(S2list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==3){
        sout$t2.surv$lower<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)-
                                           zstar*sout$t2.surv$std.err,0))
        sout$t2.surv$upper<- pmin(1,pmax(apply(S2list,2,mean,na.rm=TRUE)+
                                           zstar*sout$t2.surv$std.err,0))
      }

    }
  }
  return(sout)
}

#' @aliases scrsurvse_noJ
#' @noRd
#'
scrsurvse_noJ<- function(S1fit,S2fitkm,method,conftype,se.method,B,
                         time1,event1,time2,event2,subset=NULL,
                         copulaparam,copulaparams=NULL,bootind=NULL,copulafam,model,
                         casewt,zetas=NULL,zstar){
  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(casewt)){
    casewt<- casewt[subset]
  }


  S1list<- list()
  H1list<- list()

  n<- length(time1)
  if(conftype %in%1:2){
    if(method=="L"){
      if(se.method=="jackknife"){
        for (k in 1:n) {
          fit<- scrsurvfit.LRA(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],event2 = event2[-k],
                               copulaparam = copulaparam,copulafam = copulafam,model = model,casewt=casewt[-k])
          S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
          H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
        }
      }else{
        for (k in 1:B) {
          indd<- sort(sample(1:n,replace = TRUE))
          fit<- scrsurvfit.LRA(time1 = time1[indd],event1 = event1[indd],
                               time2 = time2[indd],event2 = event2[indd],
                               copulaparam = copulaparam,copulafam = copulafam,
                               model = model,casewt=casewt[indd])
          S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
          H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
        }
      }

      S1list<- do.call(rbind,S1list)
      H1list<- do.call(rbind,H1list)
      S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
      S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
      if(conftype==2){
        S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==1){
        S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                    zstar*S1fit$std.err,0))
        S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                    zstar*S1fit$std.err,0))
      }

    }else{

      if(se.method=="jackknife"){
        for (k in 1:n) {

          fit2<- survival::survfit( survival::Surv(time2[-k],event2[-k])~1,model = model)

          fit1<- switch (method,
                         "F" = scrsurvfit.FJC(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],event2 = event2[-k],
                                              copulaparam = copulaparam,copulafam = copulafam,
                                              model = model,casewt=casewt[-k],S2=fit2),
                         "K" = scrsurvfit.kernel(time1 = time1[-k],event1 = event1[-k],time2 = time2[-k],
                                                 event2 = event2[-k], S2=fit2)
          )
          S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
          H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)

        }
      }else{
        for (k in 1:B) {
          indd<- sort(sample(1:n,replace = TRUE))
          fit2<- survival::survfit(survival::Surv(time2[indd],event2[indd])~1,model = model)

          fit1<- switch (method,
                         "F" = scrsurvfit.FJC(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],event2 = event2[indd],
                                              copulaparam = copulaparam,copulafam = copulafam,
                                              model = model,casewt=casewt[indd],S2=fit2),
                         "K" = scrsurvfit.kernel(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],
                                                 event2 = event2[indd],S2=fit2)
          )
          S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
          H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)

        }
      }

      S1list<- do.call(rbind,S1list)
      H1list<- do.call(rbind,H1list)
      S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
      S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
      if(conftype==2){
        S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==1){
        S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                    zstar*S1fit$std.err,0))
        S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                    zstar*S1fit$std.err,0))
      }

    }
  }
  if(conftype %in%3:4){
    if(method=="L"){
      for (k in 1:B) {
        if(is.null(bootind)){
          indd<- sort(sample(1:n,replace = TRUE))
        }else{
          indd<- bootind[k,]
        }
        fit<- scrsurvfit.LRA(time1 = time1[indd],event1 = event1[indd],
                             time2 = time2[indd],event2 = event2[indd],
                             copulaparam = copulaparams[k],copulafam = copulafam,
                             model = model,casewt=casewt[indd])
        S1list[[k]]<- stepfun(fit$time,c(1,fit$surv),right = FALSE)(S1fit$time)
        H1list[[k]]<- stepfun(fit$time,c(0,fit$cumhaz),right = FALSE)(S1fit$time)
      }

      S1list<- do.call(rbind,S1list)
      H1list<- do.call(rbind,H1list)
      S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
      S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
      if(conftype==4){
        S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==3){
        S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                    zstar*S1fit$std.err,0))
        S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                    zstar*S1fit$std.err,0))
      }

    }else{

      for (k in 1:B) {
        if(is.null(bootind)){
          indd<- sort(sample(1:n,replace = TRUE))
        }else{
          indd<- bootind[k,]
        }
        fit2<- survival::survfit(survival::Surv(time2[indd],event2[indd])~1,model = model)

        fit1<- switch (method,
                       "F" = scrsurvfit.FJC(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],event2 = event2[indd],
                                            copulaparam = copulaparams[k],copulafam = copulafam,
                                            model = model,casewt=casewt[indd],S2=fit2),
                       "K" = scrsurvfit.kernel(time1 = time1[indd],event1 = event1[indd],time2 = time2[indd],
                                               event2 = event2[indd],S2=fit2)
        )
        S1list[[k]]<- stepfun(fit1$time,c(1,fit1$surv),right = FALSE)(S1fit$time)
        H1list[[k]]<- stepfun(fit1$time,c(0,fit1$cumhaz),right = FALSE)(S1fit$time)

      }


      S1list<- do.call(rbind,S1list)
      H1list<- do.call(rbind,H1list)
      S1fit$std.err<- apply(S1list,2,sd,na.rm=TRUE)
      S1fit$std.chaz<- apply(H1list,2,sd,na.rm=TRUE)
      if(conftype==4){
        S1fit$lower<- apply(S1list,2,quantile,probs = 0.05/2,na.rm=TRUE)
        S1fit$upper<- apply(S1list,2,quantile,probs = 1-0.05/2,na.rm=TRUE)
      }else if(conftype==3){
        S1fit$lower<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)-
                                    zstar*S1fit$std.err,0))
        S1fit$upper<- pmin(1,pmax(apply(S1list,2,mean,na.rm=TRUE)+
                                    zstar*S1fit$std.err,0))
      }

    }
  }

  sout<- list(t1.surv=S1fit,t2.surv=S2fitkm)
  return(sout)
}

#' @aliases scrsurvfit.FJC
#' @noRd
#'
scrsurvfit.FJC<- function(time1,event1,time2,event2,copulaparam,copulafam,model,
                          casewt,S2,subset=NULL){

  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(casewt)){
    casewt<- casewt[subset]
  }
  # ptime<- unique(as.numeric(sort(time1)))
  ptime<- as.numeric(sort(time1))
  ptime2<- unique(as.numeric(sort(time1[event1==1])))

  t2surv<- approx(c(0,S2$time),c(1,S2$surv),yright = min(S2$surv),ties = mean,
                  method = "constant",xout = ptime2,f = 0)$y


  event3<- event1+(1-event1)*event2
  t12fit <- prodlim::prodlim(survival::Surv(time1,event3)~1)

  t12surv<- predict(t12fit,times = ptime2,level.chaos='1',type="surv",mode="matrix",bytime=TRUE)


  surv1<- sapply(archmCopulaLink(copulafam,param = copulaparam, p=t12surv)-
                   archmCopulaLink(copulafam,param = copulaparam, p=t2surv),
                 FUN = archmCopulaLink_inv,copulafam = copulafam,param = copulaparam)
  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0

  F1<- quantreg::rearrange(stepfun(x= ptime2,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1

  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0


  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))


  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
               surv=surv1,cumhaz=cumhaz,type="right")


  return(out)

}


#' @aliases scrsurvfit.LRA
#' @noRd
#'
scrsurvfit.LRA<- function(time1,event1,time2,event2,copulaparam,copulafam,
                         model,casewt,subset=NULL){

  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(casewt)){
    casewt<- casewt[subset]
  }

  ptime<- as.numeric(sort(time1))
  ptime2<- unique(as.numeric(sort(time1[event1==1])))

  event3<- event1+(1-event1)*event2
  t12fit <- prodlim::prodlim(survival::Surv(time1,event3)~1)

  t12surv<- predict(t12fit,times = ptime,level.chaos='1',type="surv",mode="matrix",bytime=TRUE)

  surv1<- sapply(ptime2,function(t){
    pos<- max(which(ptime==t))
    archmCopulaLink_inv(copulafam,param = copulaparam,
                        y=sum(archmCopulaLink(copulafam,param = copulaparam, p=t12surv[1:pos])-
                                archmCopulaLink(copulafam,param = copulaparam, p=c(1,t12surv)[1:pos]),
                              na.rm = TRUE))
  })

  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0

  ptime<- unique(ptime)

  F1<- quantreg::rearrange(stepfun(x= ptime2,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1

  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0



  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))


  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
               surv=surv1,cumhaz=cumhaz,type="right")



  return(out)

}


#' @aliases scrsurvfit.JFKC
#' @noRd
scrsurvfit.JFKC<-function(time1,event1,time2,event2,copulaparam,copulafam,
                          model,casewt,surv2km=FALSE,subset=NULL){

  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  if(!is.null(subset)&!is.null(casewt)){
    casewt<- casewt[subset]
  }
  tp<- sort(unique(c(time1[event1==1],time2[event2==1])))
  # tp<- sort(unique(c(time1,time2)))
  term1 <- sapply(tp,function(t) sum((time1> t)*casewt,na.rm=TRUE))/sum(casewt)
  if(!isTRUE(surv2km)){
    term2 <- sapply(tp,function(t) sum((time2> t)*casewt,na.rm=TRUE))/sum(casewt)
  }

  S1<- survival::survfit(survival::Surv(time1,event1)~1,weights = casewt)
  F1bar.ini<- stepfun(x = c(0,S1$time),y= c(1,S1$surv,0),right = FALSE,ties = mean)(tp)

  S2<- survival::survfit(survival::Surv(time2,event2)~1,weights = casewt)
  F2bar.ini<- stepfun(x = c(0,S2$time),y= c(1,S2$surv,0),right = FALSE,ties = mean)(tp)

  it<- 1
  mit<- 200
  test0<- 1000
  convergence<- 0
  repeat{
    F1bar.ini.x<- approx(x=tp,y=F1bar.ini,xout = time1,method = "constant",f = 0,ties = mean)$y
    F2bar.ini.x<- approx(x=tp,y=F2bar.ini,xout = time2,method = "constant",f = 0,ties = mean)$y

    dc10<- sapply(1:length(tp),function(k)
      sum(Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini[k],p2 = F2bar.ini.x)/
            Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini.x)*
            (1-event1)*(1-event2)*casewt*(time1<= tp[k]),na.rm = TRUE)/sum(casewt))

    dc12<- sapply(1:length(tp),function(k)
      sum(dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini[k],p2 = F2bar.ini.x,mode = "2")/
            dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini.x,mode = "2")*
            (1-event1)*event2*casewt*(time1<= tp[k]),na.rm = TRUE)/sum(casewt))
    dc10[is.na(dc10)]<- 0
    dc10[is.infinite(dc10)]<- 0

    dc12[is.na(dc12)]<- 0
    dc12[is.infinite(dc12)]<- 0
    F1bar.est<- term1 + dc10+ dc12

    if(!isTRUE(surv2km)){
      dc20<- sapply(1:length(tp),function(k)
        sum(Copulafn(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini[k])/
              Copulafn(copulafam = copulafam,param = copulaparam, p1 = F1bar.ini.x,p2 = F2bar.ini.x)*
              (1-event1)*(1-event2)*casewt*(time2<= tp[k]),na.rm = TRUE)/sum(casewt))

      dc21<- sapply(1:length(tp),function(k)
        sum(dev_Copula(copulafam = copulafam,param = copulaparam,p1 = F1bar.ini.x,p2 = F2bar.ini[k],mode = "1")/
              dev_Copula(copulafam = copulafam,param = copulaparam, p1 = F1bar.ini.x,p2 = F2bar.ini.x,mode = "1")*
              event1*(1-event2)*casewt*(time2<= tp[k]),na.rm = TRUE)/sum(casewt))
      dc20[is.na(dc20)]<- 0
      dc20[is.infinite(dc20)]<- 0

      dc21[is.na(dc21)]<- 0
      dc21[is.infinite(dc21)]<- 0
      F2bar.est<- term2 + dc20 + dc21
      F2bar.est<- pmin(pmax(0,F2bar.est),1)
    }else{
      F2bar.est<-F2bar.ini
    }

    F1bar.est<- pmin(pmax(0,F1bar.est),1)


    test<- max(max(abs(F1bar.est-F1bar.ini)),max(abs(F2bar.est-F2bar.ini)))
    # if(it%%20==1){cat(it,test,".\n")}
    if(test< 1e-4|it>mit){
      # cat(it,test,".\n")
      convergence=0
      break
    }else if(test>test0){
      # cat(it,test,".\n")
      F1bar.est <- F1bar.ini
      F2bar.est <- F2bar.ini
      if(test>0.05){convergence=1}
      break
    }else{
      it<-it+1
      F1bar.ini<- F1bar.est
      F2bar.ini<- F2bar.est
      test0<- test
    }

  }



  F1<- 1- F1bar.est
  F2<- 1- F2bar.est

  napos<- union(which(is.na(F1)),which(is.na(F2)))
  if(length(napos)>0){
    tp<- tp[-napos]
    F1<- F1[-napos]
    F2<- F2[-napos]

  }

  # ptime1<-  as.numeric(sort(unique(time1)))
  # ptime2<-  as.numeric(sort(unique(time2)))
  ptime1<-  as.numeric(sort(time1))
  ptime2<-  as.numeric(sort(time2))
  F1<- quantreg::rearrange(stepfun(x= tp,y = c(0,F1),right = FALSE,ties = mean))(ptime1)
  F2<- quantreg::rearrange(stepfun(x= tp,y = c(0,F2),right = FALSE,ties = mean))(ptime2)


  surv1<- 1-F1
  surv2<- 1-F2

  n.risk1<- as.numeric(sapply(ptime1,function(tt)sum(time1>=tt)))
  n.event1<- as.numeric(sapply(ptime1,function(tt)sum(event1*(time1==tt))))
  n.censor1<- as.numeric(sapply(ptime1,function(tt)sum((1-event1)*(time1==tt))))

  cumhaz1<-  -log(surv1)
  cumhaz1[cumhaz1<0 ]<- 0



  out1<- list(time=ptime1,n.risk=n.risk1,n.event=n.event1,n.censor=n.censor1,
              surv=surv1,cumhaz=cumhaz1,type="right")
  if(!isTRUE(surv2km)){
    n.risk2<- as.numeric(sapply(ptime2,function(tt)sum(time2>=tt)))
    n.event2<- as.numeric(sapply(ptime2,function(tt)sum(event2*(time2==tt))))
    n.censor2<- as.numeric(sapply(ptime2,function(tt)sum((1-event2)*(time2==tt))))

    cumhaz2<-  -log(surv2)
    cumhaz2[cumhaz2<0 ]<- 0
    out2<- list(time=ptime2,n.risk=n.risk2,n.event=n.event2,n.censor=n.censor2,
                surv=surv2,cumhaz=cumhaz2,type="right")
  }else{
    out2<- S2
  }

  return(list(s1=out1,s2=out2))




}


#' @aliases scrsurvfit.kernel
#' @noRd
scrsurvfit.kernel<- function(time1,event1,time2,event2,S2,subset=NULL){

  # A modification of codes from Nevo D, Gorfine M. (2022)
  # Causal inference for semi-competing risks data.
  # Biostatistics, 14, 23(4):1115-1132.
  if(!is.null(subset)){
    time1<- time1[subset]
    event1<- event1[subset]
    time2<- time2[subset]
    event2<- event2[subset]

  }
  # if(!is.null(subset)&!is.null(casewt)){
  #   casewt<- casewt[subset]
  # }

  # ptime<- unique(as.numeric(sort(time1)))
  ptime<- as.numeric(sort(time1))
  n.times<- length(ptime)

  dS2 <-  diff(c(1, S2$surv))
  dS2.times <-  S2$time
  dS2.times<- dS2.times[dS2!= 0]
  dS2 <- dS2[dS2!= 0]
  n.times2 <- length(dS2.times)

  T1dead <- time1[event2==1]
  event1dead <- event1[event2==1]
  T2dead <- time2[event2==1]

  # Calculate S_{1|A=0, T_2=t}
  S1_2 <- prodlim::prodlim(survival::Surv(T1dead, event1dead) ~ T2dead)
  newdata <- data.frame(T2dead = dS2.times)
  S1_2 <- predict(S1_2, newdata = newdata, times =  dS2.times)
  S1_2 <- do.call(rbind, S1_2)


  # Change NA after no more events into the last estimated value
  S1_2 <- t(apply(S1_2, 1, RepNAmin))

  S1.ptime <-  vector(length = n.times)
  S1st <- matrix(nrow = n.times, ncol = n.times2)

  S1AtT2s <- function(my.t, my.s) {
    if (my.s < my.t){      my.t <- my.s}
    s.place <- findInterval(my.s, dS2.times)
    t.place <- findInterval(my.t, dS2.times)
    if (t.place==0 | s.place==0) {
      S.res <- 1
    } else {
      t1.curve.s <- S1_2[s.place, ]
      if (t.place==length(dS2.times)) {
        S.res <-  t1.curve.s[length(dS2.times)]
      } else {
        S.res <- t1.curve.s[t.place]
      }}

    return(S.res)
  }
  for (j in 1:n.times) {
    t.now <- ptime[j]
    pos.t.now <- findInterval(t.now, dS2.times)

    if(pos.t.now==0) {
      S1.ptime[j] <- 1
    } else {
      S1st[j, ] <-  sapply(dS2.times, S1AtT2s, my.t = t.now)
      S1.ptime[j] <- -sum(S1st[j, ] * dS2)
    }



  }


  surv1<- S1.ptime
  surv1[surv1>1]<- 1
  surv1[surv1<0]<- 0
  F1<- quantreg::rearrange(stepfun(x= ptime,y = c(0,1-surv1),right = FALSE,ties = mean))(ptime)
  surv1<- 1-F1

  cumhaz<-  -log(surv1)
  cumhaz[cumhaz<0 ]<- 0


  n.risk<- as.numeric(sapply(ptime,function(tt)sum(time1>=tt)))
  n.event<- as.numeric(sapply(ptime,function(tt)sum(event1*(time1==tt))))
  n.censor<- as.numeric(sapply(ptime,function(tt)sum((1-event1)*(time1==tt))))


  out<- list(time=ptime,n.risk=n.risk,n.event=n.event,n.censor=n.censor,
             surv=surv1,cumhaz=cumhaz,type="right")

  return(out)
}


#' @aliases RepNAmin
#' @noRd
RepNAmin <- function(x){
  min.x <- min(x, na.rm = T)
  x[is.na(x)] <- min.x
  return(x)
}

#' @aliases survtransform
#' @noRd
surv_transform<- function(survlist,inverse=FALSE,conf.int,t1data,t2data){

  if(!isTRUE(inverse)){
    size.strata<- t1data$size.strata
    ## when ns>1 preform transformation for surv list
    ns<- length(survlist)
    newlist<- list()

    surv1<- lapply(survlist,function(x)x$t1.surv)

    newlist$t1.surv$time<- do.call("c",lapply(surv1,function(x)x$time))
    newlist$t1.surv$n.risk<- do.call("c",lapply(surv1,function(x)x$n.risk))
    newlist$t1.surv$n.event<- do.call("c",lapply(surv1,function(x)x$n.event))
    newlist$t1.surv$n.censor<- do.call("c",lapply(surv1,function(x)x$n.censor))
    newlist$t1.surv$surv<- do.call("c",lapply(surv1,function(x)x$surv))
    newlist$t1.surv$cumhaz<- do.call("c",lapply(surv1,function(x)x$cumhaz))

    newlist$t1.surv$first.strata<- cumsum(c(1,size.strata))[1:length(size.strata)]
    newlist$t1.surv$size.strata<- size.strata
    newlist$t1.surv$X<- t1data$stratas
    # newlist$t1.surv$maxtime<- max(newlist$t1.surv$time)

    newlist$t1.surv$discrete.predictors<- t1data$discrete.predictors
    newlist$t1.surv$continuous.predictors<- t1data$continuous.predictors
    newlist$t1.surv$xlevels<- t1data$xlevels


    surv2<- lapply(survlist,function(x)x$t2.surv)

    newlist$t2.surv$time<- do.call("c",lapply(surv2,function(x)x$time))
    newlist$t2.surv$n.risk<- do.call("c",lapply(surv2,function(x)x$n.risk))
    newlist$t2.surv$n.event<- do.call("c",lapply(surv2,function(x)x$n.event))
    newlist$t2.surv$n.censor<- do.call("c",lapply(surv2,function(x)x$n.censor))
    newlist$t2.surv$surv<- do.call("c",lapply(surv2,function(x)x$surv))
    newlist$t2.surv$cumhaz<- do.call("c",lapply(surv2,function(x)x$cumhaz))

    newlist$t2.surv$first.strata<- cumsum(c(1,size.strata))[1:length(size.strata)]
    newlist$t2.surv$size.strata<- size.strata

    # newlist$t2.surv$maxtime<- max(newlist$t2.surv$time)
    newlist$t2.surv$X<- t2data$stratas

    newlist$t2.surv$discrete.predictors<- t2data$discrete.predictors
    newlist$t2.surv$continuous.predictors<- t2data$continuous.predictors
    newlist$t2.surv$xlevels<- t2data$xlevels



    if(is.numeric(conf.int)){

      newlist$t1.surv$conf.int<- newlist$t2.surv$conf.int<- conf.int

      newlist$t1.surv$std.err<- do.call("c",lapply(surv1,function(x)x$std.err))
      newlist$t1.surv$std.chaz<- do.call("c",lapply(surv1,function(x)x$std.chaz))
      newlist$t1.surv$lower<- do.call("c",lapply(surv1,function(x)x$lower))
      newlist$t1.surv$upper<- do.call("c",lapply(surv1,function(x)x$upper))


      newlist$t2.surv$std.err<- do.call("c",lapply(surv2,function(x)x$std.err))
      newlist$t2.surv$std.chaz<- do.call("c",lapply(surv2,function(x)x$std.chaz))
      newlist$t2.surv$lower<- do.call("c",lapply(surv2,function(x)x$lower))
      newlist$t2.surv$upper<- do.call("c",lapply(surv2,function(x)x$upper))

    }
  }



  return(newlist)

}


