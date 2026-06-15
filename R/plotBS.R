#' Brier Scores for terminal survival prediction
#' @aliases dyBS
#' @usage dyBS(surv, S_D, T2, event2, times=knots(S_D), reference = TRUE,
#'   tu=NULL, int.method=c("none","fmm","natural","periodic","monoH.FC","hyman"))
#' @param surv a list containing \code{surv} (survival prediction output from object of class "mscr") and \code{t1} (available observations for multiple nonterminal event times).
#' @param S_D (estimated or specified) survival function for the terminal event time. See Examples. 
#' @param T2 A numeric vector of observed terminal event times. Its length must
#'   equal the number of rows of \code{event2}.
#' @param event2 A numeric or logical vector of censoring indicators for the
#'   terminal event time, where 1 indicates the terminal event is observed and
#'   0 indicates censoring.
#' @param times evaluation times for the Brier score curve
#' @param reference If TRUE, the returned object is a list containing the relative Brier Score compared to the predicted survival curve yielded from the benchmark Kaplan Meier estimator.
#' @param tu A numeric value specifying the maximum time point for the Brier score curve. If \code{tu = NULL}, no maximum time point is used.
#' @param int.method	specifies the type of spline to be used.
#'   Possible values are "none","fmm", "natural", "periodic", "monoH.FC" and "hyman". 
#'   See also \code{\link[stats]{splinefun}}.
#' @return A list containing the following components: BS=data.frame(time=times,BS=BS), IBS=IBS
#' \describe{
#'     \item{BS}{Estimated Brier Score Curve.}
#'     \item{IBS}{Integrated Brier Score over the range of evaluation times.}
#'   }
#' @export
#'
#' @seealso \code{\link{mscr}}, \code{\link{predict.mscr}}
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' data<- simSCRmul(n=100,copulafam = "frank",K=3,
#'                  tau.alpha = .5,tau.theta=.5)
#'
#' fit<- mscr(mT1=data$T1,T2=data$T2,mevent1=data$event1,event2=data$event2,
#'            copulafam = "frank")
#' mT1.te<- data$T1
#' mT1.te[data$event1==0]<- NA
#'
#' ####  dynamic prediction
#' surv<- predict(fit, t1obs=mT1.te,type="surv")
#'
#' ## naive Pr(D>t|D> t_m)
#' surv0<- predict(fit, t1obs=mT1.te,type="surv",cause=0)
#'
#' ## Pr(D>t|T_1=t_1,D> t_1)
#' surv1<- predict(fit, t1obs=mT1.te,type="surv",cause=1,maxobs = FALSE)
#'
#' ## Pr(D>t|T_1=t_1,D> t_m)
#' surv1m<- predict(fit, t1obs=mT1.te,type="surv",cause=1,maxobs = TRUE)
#'
#' ### Brier Scores
#' BS<- dyBS(surv,S_D = fit$S_D,T2=data$T2,event2=data$event2,int.method = "natural")
#' BS0<- dyBS(surv0,S_D = fit$S_D,T2=data$T2,event2=data$event2,int.method = "natural")
#' BS1<- dyBS(surv1,S_D = fit$S_D,T2=data$T2,event2=data$event2,int.method = "natural")
#' BS1m<- dyBS(surv1m,S_D = fit$S_D,T2=data$T2,event2=data$event2,int.method = "natural")
#' (ibs<-c(BS$IBS,BS0$IBS,BS1$IBS,BS1m$IBS)) ## integrated BS
#'
#' ### plot Brier Scores
#' plot(BS$BS,ylim = c(0,0.7),type="l",xlab="time",ylab="Brier Score",col=1)
#' lines(BS0$BS,col=2)
#' lines(BS1$BS,col=3)
#' lines(BS1m$BS,col=4)
#' legend("topright",col = 1:4,lty = 1,legend = c("DP","P0","P1","P1m"))
#' }


dyBS<- function(surv,S_D,T2,event2,times=knots(S_D),
                reference=TRUE,tu=NULL,
                int.method=c("none","fmm", "natural", "periodic", "monoH.FC","hyman")){

  int.method<- match.arg(int.method)
  times<- unique(c(0,times))

  if(!is.null(tu)){times<- times[times<= tu]}
  if(isTRUE(reference)){
    BS00<- baseBS(S_D=S_D,T2=T2,event2=event2,times=times)$BS$BS
  }


  cens.fit = prodlim::prodlim(Surv(T2,event2)~1,
                              data = data.frame(T2 = T2,event2 =event2),reverse=TRUE)
  ipcw = predict(cens.fit,times=times,level.chaos=1,mode="matrix",type="surv")
  ipcw.obs = prodlim::predictSurvIndividual(cens.fit,lag=1)
  ipcw[ipcw==0]=min(ipcw[ipcw!=0],2e-10)/2
  ipcw.obs[ipcw.obs==0]=min(ipcw.obs[ipcw.obs!=0],2e-10)/2

  SurvMatrix<- lapply(surv$surv,function(x)
    if(sum(!is.na(x$surv))>0){
      stats::approx(x=x$times,y=x$surv,xout = times,method = "constant",ties = mean,
             yleft = 1,yright = min(x$surv,na.rm = TRUE),f = 0,na.rm = TRUE)$y
    }else{
      rep(NA,length(times))
    })


  SurvMatrix2<- do.call(rbind,SurvMatrix)
  BS=sapply(times, function(t){
    mean(as.numeric(T2<=t)*event2*(SurvMatrix2[,times==t]^2)/ipcw.obs,na.rm =TRUE)+
      mean(as.numeric(T2>t)*(1-SurvMatrix2[,times==t])^2/ipcw[times==t],na.rm =TRUE)
  })
  if(int.method=="none"){
    IBS=1/diff(range(times))*sum(BS*diff(c(times,max(times))),na.rm =TRUE)
  }else{
    fn<- stats::splinefun(times[!is.na(BS)],BS[!is.na(BS)],method = int.method)
    IBS<- pracma::quadinf(f = fn,xa = min(times),xb = max(times))$Q
  }


  if(isTRUE(reference)){
    RBS<- 1- BS/BS00
    RBS[is.na(RBS)|is.infinite(RBS)]<- 0
    IRBS=1/diff(range(times))*sum(RBS*diff(c(times,max(times))),na.rm =TRUE)

    return(list(BS=data.frame(time=times,BS=BS),
                RBS=data.frame(time=times,RBS=RBS),
                IBS=IBS,IRBS=IRBS))

  }else{
    return(list(BS=data.frame(time=times,BS=BS),IBS=IBS))
  }


}


#' @noRd
baseBS<- function(S_D,T2,event2,times=knots(S_D)){
  cens.fit = prodlim::prodlim(Surv(T2,event2)~1,
                              data = data.frame(T2 = T2,event2 =event2),reverse=TRUE)
  ipcw = predict(cens.fit,times=times,level.chaos=1,mode="matrix",type="surv")
  ipcw.obs = prodlim::predictSurvIndividual(cens.fit,lag=1)
  ipcw[ipcw==0]=min(ipcw[ipcw!=0],2e-10)/2
  ipcw.obs[ipcw.obs==0]=min(ipcw.obs[ipcw.obs!=0],2e-10)/2

  SurvMatrix<- S_D(times)

  SurvMatrix2<- matrix(SurvMatrix,byrow = TRUE,nrow = length(T2),ncol = length(times))
  BS=sapply(times, function(t){
    mean(as.numeric(T2<=t)*event2*(SurvMatrix2[,times==t]^2)/ipcw.obs,na.rm =TRUE)+
      mean(as.numeric(T2>t)*(1-SurvMatrix2[,times==t])^2/ipcw[times==t],na.rm =TRUE)
  })


  IBS=1/diff(range(times))*sum(BS*diff(c(times,max(times))),na.rm =TRUE)
  return(list(BS=data.frame(time=times,BS=BS),IBS=IBS))

}
