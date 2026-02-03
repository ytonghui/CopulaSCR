#################################################################
##   CopulaSCR R package by Tonghui Yu Copyright (C) 2024
##
## The below functions are used for survival prediction for Semi-competing
## risks data with multiple intermediate event times
#################################################################
#' Survival prediction
#'
#' Estimation of  joint survival function, cross-hazard ratio function
#' and conditional survival function under semi-competing risks data with
#' single intermediate event.
#'
#' @name  predict.scrsurv
#' @aliases predict.scrsurv
#' @usage predict(object, newdata, type = c("msurv","jsurv","csurv","chr"),
#'    t1, t2, s2, t1equal = FALSE, bystrata = TRUE,...)
#' @param object an optional object of class  "scrsurv".
#' @param newdata data frame in which to interpret the variables
#'   occurring in the formula.
#' @param type specifies the type of function to be estimated.
#' @param t1 numeric vector
#' @param t2 numeric vector
#' @param s2 numeric vector
#' @param t1equal If TRUE,..
#' @param bystrata If TRUE,...
#' @param ... other arguments
#' @return A list with components: t1, t2, and one of jsurv, msurv, csurv, and chr.
#' @export
#' @seealso \code{\link{scrassonp}}, \code{\link{scrsurv}}, \code{\link{predictscr}}
#' @keywords survival
#' @examples
#' library(CopulaSCR)
#' set.seed(12345)
#' simdata<- simSCR(n = 50,tau = 0.5, copulafam = "frank")
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
#'                     t2.formula=Surv(T2, event2) ~ 1,
#'                     data = simdata,copulafam="frank",
#'                     a=quantile(simdata$T1,0.9),
#'                     b=quantile(simdata$T2,0.9),B=10,seed = 12345,
#'                     se = TRUE,se.method = "resampling")
#' summary(fitasso)
#' fit<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=fitasso$B,
#' se.method = fitasso$se.method,conf.int=TRUE,conftype =3)
#'
#' t1<- t2<- seq(0,3.5,0.1)
#' ### marginal survival probabilities
#' ms<- predict(fit,t1=t1, type="msurv")
#'
#' ### joint survival probabilities
#' js<- predict(fit,t1=t1, t2=t2,type="jsurv")
#'
#' ## stratification analysis
#' set.seed(12345)
#' simdata2<- simSCRtr(n = 100,K=3, tau = c(0.3,0.5,0.6), copulafam = "frank",
#'                     params=list(marginsDist = rep("exp",2), rate1=c(0.5,1,1.2),rate2=1))
#' fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ tr,
#'                     t2.formula=Surv(T2, event2) ~ tr,
#'                     data=simdata2, copulafam="frank",a=quantile(simdata2$T1,0.9),
#'                     b=quantile(simdata2$T2,0.9),B=0)
#' fittr<- scrsurv(fit = fitasso, method= "JFKC",surv2km=FALSE,B=0)
#' t1<- s2<- c(0.2,0.5,1)
#' t2<-  seq(0,3.5,0.1)
#' ms2<- predict(fittr,t1=t1, type="msurv",newdata=simdata2[1:3,])
#' js2<- predict(fittr,t1=t1, t2=t2, type="jsurv",newdata=simdata2[1:3,])
#' cs2<- predict(fittr,t1=t1, t2=t2,s2=s2, type="csurv",newdata=simdata2[1:3,])
#' chr2<- predict(fittr,t1=t1, t2=t2, type="chr",newdata=simdata2[1:3,])
#'

"predict.scrsurv"<- function(object, newdata,type = c("msurv","jsurv","csurv","chr"),
                           t1, t2, s2, t1equal = FALSE,bystrata = TRUE,...){
  Call <- match.call()
  type <- match.arg(type)
  indx <- match(c('object','t1', 't2' ,'s2','newdata'),names(Call), nomatch=0)
  if(type!="msurv"){
    if (indx[2]==0) stop("a t1 argument is required")
    if (indx[3]==0) stop("a t2 argument is required")
  }
  if (type=="csurv"&indx[4]==0) stop("a s2 argument is required")

  if(object$t1.surv$covariate.type!=1){
    if (indx[5]==0) stop("a newdata argument is required")
    if (!(is.matrix(newdata) | is.data.frame(newdata))) {
      stop("data object must be a matrix or a data.frame")
    }
  }

  if(missing(newdata))newdata<- NULL
  if(type=="msurv"){
    if (indx[2]==0&indx[3]==0) stop("either of t1 and t2 arguments is required")
    out<- list()

    if (object$t1.surv$covariate.type==1){

      if(indx[2]!=0){
        if(is.null(t1))stop("a t1 argument is required")
        psurv1<- stepfun(x= object$t1.surv$time,y=c(1,object$t1.surv$surv),
                         right = FALSE,ties = max)
        out$t1<- t1
        out$msurv1<- psurv1(t1)
      }
      if(indx[3]!=0){
        if(is.null(t2))stop("a t2 argument is required")
        psurv2<- stepfun(x= object$t2.surv$time,y=c(1,object$t2.surv$surv),
                         right = FALSE,ties = max)
        out$t2<- t2
        out$msurv2<- psurv2(t2)
      }

    }  else{
      strata.vars <- object$t1.surv$discrete.predictors
      # NN.vars <- object$t1.surv$continuous.predictors
      X.formula <- update(formula(object$t1.formula),NULL~.)
      if (!all(match(all.vars(X.formula),names(newdata),nomatch=FALSE)))
        stop("Arg newdata does not contain all the covariates used for fitting. \n\nfitted variables: ", paste(all.vars(X.formula),collapse=", "),"\nnewdata contains:",ifelse(length(names(newdata))==0," nothing",names(newdata)))
      requested.X <- newdata[,all.vars(X.formula),drop=FALSE]


      fit.strata <- interaction(object$stratas[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      requested.strata <- interaction(requested.X[,strata.vars,drop=FALSE],sep=":",drop=TRUE)
      fit.levels <- as.character(unique(fit.strata))
      NS <- length(fit.levels)
      fit.strata <- factor(fit.strata,levels=levels(fit.strata),labels=1:NS)
      requested.strata <- factor(requested.strata,levels=fit.levels,labels=1:NS)
      strata<- rep(fit.strata,times=object$t1.surv$size.strata)

      if(indx[2]!=0){
        if(is.null(t1))stop("a t1 argument is required")
        survmat1<- data.frame(time=object$t1.surv$time,surv=object$t1.surv$surv,strata=strata)


        survfn1<- lapply(fit.strata,function(k){
          sv<- survmat1[survmat1[,3]==k,]
          stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
        })
        psurv1<- survfn1[requested.strata]
        psurv1<- lapply(psurv1,function(x)x(t1))
        if (bystrata == TRUE){
          names(psurv1) <- paste(strata.vars,requested.X[,1],sep="=")
        }else{
          psurv1<- do.call(cbind,psurv1)
          psurv1<- lapply(seq(1:length(t1)),function(k)as.numeric(psurv1[k,]))
          names(psurv1) <- paste("t1",t1,sep="=")
        }

        out$t1<- t1
        out$msurv1<- psurv1
      }

      if(indx[3]!=0){
        if(is.null(t2))stop("a t2 argument is required")
        survmat2<- data.frame(time=object$t2.surv$time,surv=object$t2.surv$surv,strata=strata)
        survfn2<- lapply(fit.strata,function(k){
          sv<- survmat2[survmat2[,3]==k,]
          stepfun(x= sv[,1],y=c(1,sv[,2]),right = FALSE,ties = max)
        })
        psurv2<- survfn2[requested.strata]
        psurv2<- lapply(psurv2,function(x)x(t2))
        if (bystrata == TRUE){
          names(psurv2) <- paste(strata.vars,requested.X[,1],sep="=")
        }else{
          psurv2<- do.call(cbind,psurv2)
          psurv2<- lapply(seq(1:length(t2)),function(k)as.numeric(psurv2[k,]))
          names(psurv2) <- paste("t2",t2,sep="=")
        }
        out$t2<- t2
        out$msurv2<- psurv2
      }
    }
  }

  if(type=="jsurv"){
    out<- predjsurv(fit=object, t1=t1, t2=t2, method = "sp",newdata =newdata)
  }
  if(type=="csurv"){
    out<- predcsurv(fit=object, t1=t1, t2=t2, s2=s2, t1equal = t1equal,
                    method = "sp",newdata=newdata)
  }
  if(type=="chr"){
    out<- predchr(fit=object, t1=t1, t2=t2, method = "sp",newdata =newdata)
  }
  return(out)
}


#' Dynamic terminal survival prediction
#'
#' Fitting semi-competing risks data with multiple intermediate event times using
#' a copula-based model. Survival prediction for terminal event based on the observed
#' intermediate event times.
#'
#' @name predict.mscr
#' @aliases predict.mscr
## #' @usage predictmscr(fit, t1obs, times, type = c("surv", "rrms", "rmst", "qrl",
## #'   "qst"), tu, tau=0.5, nsim = 1000, smooth = FALSE, jumpsize = FALSE, exact = FALSE)
## #' @usage predbsmscr(fit, t1obs, times, type=c("surv","rrms", "rmst", "qrl",
## #'   "qst"), tu, tau = 0.5, cause = 0, maxobs=TRUE)
## #'
#' @usage predict(fit, t1obs,times= knots(fit$S_D), cause = NULL,
#'      type=c("msurv","surv", "rrms", "rmst", "qrl","qst"),
#'      tu=max(knots(fit$S_D)),tau=0.5, nsim=1000, maxobs=TRUE,
#'      smooth = FALSE, jumpsize = FALSE, exact = TRUE,...)
#'
#' @param fit An object of class "mscr".
#' @param t1obs The observed value for intermediate event times.
#' @param times Evaluated time.
#' @param type type=c("surv", "rrms", "rmst", "qrl","qst")
#' @param tu A optional value for the maximum follow-up time.
#' @param tau Quantile level. \code{tau} is numeric value greater than 0 and smaller than 1.
#' @param nsim  Monte Carlo sample size for differential/integral calculus.
#' @param smooth Whether or not use smoothing.
#' @param jumpsize Whether or not use smoothing.
#' @param exact Whether or not use smoothing.
#' @param cause Integer. The survival prediction based on the 1st intermediate event if \code{cause=1}
#' @param maxobs TRUE or FALSE
#' @param ... other arguments
#' @return a
#' @keywords survival
#' @seealso \code{\link{mscr}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' data<- simSCRmul(n=100,copulafam = "frank",K=3,
#'                  tau.alpha = .5,tau.theta=.5)
#'
#' fit<- mscr(mT1=data$T1,T2=data$T2,mevent1=data$event1,event2=data$event2,
#'            copulafam = "frank",ncore = 1)
#' mT1<- data$T1
#' mT1[data$event1==0]<- NA
#' mT1.te<- mT1[1:4,]
#' msurv<- predict(fit, t1obs=mT1.te,type="msurv")
#'
#' ####  dynamic prediction
#' surv<- predict(fit, t1obs=mT1.te,type="surv")
#' rmst<- predict(fit, t1obs=mT1.te,type="rmst")
#' rmst<- sapply(rmst$rmst, function(x)x$rmst)
#'
#' ## naive Pr(D>t|D> t_m)
#' surv0<- predict(fit, t1obs=mT1.te,type="surv",cause=0)
#' rmst0<- predict(fit, t1obs=mT1.te,type="rmst",cause=0)
#' rmst0<- sapply(rmst0$rmst, function(x)x$rmst)
#'
#' ## Pr(D>t|T_1=t_1,D> t_1)
#' surv1<- predict(fit, t1obs=mT1.te,type="surv",cause=1,maxobs = FALSE)
#' rmst1<- predict(fit, t1obs=mT1.te,type="rmst",cause=1,maxobs = FALSE)
#' rmst1<- sapply(rmst1$rmst, function(x)x$rmst)
#'
#' ## Pr(D>t|T_1=t_1,D> t_m)
#' surv1m<- predict(fit, t1obs=mT1.te,type="surv",cause=1,maxobs = TRUE)
#' rmst1m<- predict(fit, t1obs=mT1.te,type="rmst",cause=1,maxobs = TRUE)
#' rmst1m<- sapply(rmst1m$rmst, function(x)x$rmst)
#'
#' ## compare with the true death time
#' T2<- data$T2
#' T2[data$event2==0]<- NA
#' T2<- T2[1:4]
#' ## mean square error
#' apply(abs(cbind(T2-rmst,T2-rmst0,T2-rmst1,rmst1m)^2),2,mean)
#'}
#'
#'
"predict.mscr"<- function(object, t1obs,times= knots(object$S_D),cause = NULL,
                          type=c("msurv","surv", "rrms", "rmst", "qrl","qst"),
                          tu=max(knots(object$S_D)),tau=0.5,nsim=1000,maxobs=TRUE,
                          smooth= FALSE,jumpsize=FALSE,exact=TRUE,...){
  type<- match.arg(type)
  if(type!="msurv"){
    if(is.null(cause)){
      predictMSCR(object, t1obs,times,type,tu,tau,nsim,smooth,jumpsize,exact)
    }else{
      predbsMSCR(object, t1obs,times,type,tu,tau,cause,maxobs)
    }
  }else{
    ms<- lapply(1:length(object$mar.fits),function(k)
      predict(object=object$mar.fits[[k]],t1=t1obs[,k], type="msurv"))
    names(ms)<- names(fit$mar.fits)
    return(ms)
  }

}

