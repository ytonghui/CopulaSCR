#################################################################
##   CopulaSCR R package by Tonghui Yu Copyright (C) 2024
##
## The below functions are used for association and marginal analysis for Semi-competing
## risks data with multiple intermediate event times
#################################################################
#' @title Association analysis for Semi-competing risks data with multiple intermediate event times
#' @aliases mscr
#'
#' @usage mscr(mT1, T2, mevent1, event2, solver = c("optimize", "BBspg"),
#'   copulafam = c("frank", "clayton", "joe", "gumbel", "amh"),
#'   nsim = 500, ncore = 4, tol = 0.01, a = 0, b = 0,
#'   positive=TRUE, msurv.method = "JFKC", lower = 0.02, upper = 0.96)
#'
#' @param mT1 Observed intermediate event times in a matrix or data frame form.
#' The columns correspond to multiple events, and the column names are the names of these events.
#' @param T2 Observed terminal event time in a vector form.
#' Its length shall equal the number of rows of \code{mT1}.
#' @param mevent1 Censoring indicators for intermediate event times in a matrix or
#' data frame form. Elements in \code{mT1} and \code{mevent1} must correspond to each other.
#' @param event2 Censoring indicators for terminal event time (1= occurred, 0= censored).
#' @param solver a character string specifying the optimization solver to be used for minumum search.
#' @param copulafam A character string specifying the family of an Archimedean copula.
#' @param nsim Monte Carlo sample size for differential/integral calculus.
#' @param ncore The number of cores for parallel computing.
#' @param tol The desired accuracy.
#' @param a positive constant to dampen weight function.
#' If \code{a = NULL} choose \code{quantile(T1,0.95)} for \eqn{a}
#' Default is \code{a=0}.
#' @param b positive constant dampen dampen weight function.
#' If \code{b = NULL} choose \code{quantile(T2,0.95)} for \eqn{b}.
#' Default is \code{b=0}. When both \code{a=0} and \code{b=0}, the the weight function \eqn{w = 1},
#' namely, \code{equalweight=TRUE}. (See also \code{scrassonp()})
#' @param positive whether or not put positive constraint on the marginal Kendall' tau parameter.
#' (See also \code{scrassonp()})
#' @param msurv.method a character string specifying the method for computing marginal survival
#'   curves of intermediate events. Currently supported methods include "JFKC", "LRA", "FJC", and "Ker". (See also \code{scrassonp()})
#' @param lower The lower bound of Kendall's tau correspond to copula parameter \eqn{\alpha}
#' in joint distribution of intermediate event times (see more in Details).
#' @param upper The upper bound of Kendall's tau correspond to copula parameter
#' \eqn{\alpha}.
#'
#' @return  An object of class "\code{mscr}" representing the fit.
#' The \code{mscr} object is a list containing at least the following components:\cr
#' \describe{
#' \item{alpha}{Estimated copula parameter.}
#' \item{logLik.alpha}{Value of pseudo-likelihood function evaluated at alpha.}
#' \item{alpha.int}{Initial end-points of the interval to be searched for alpha.}
#' \item{call}{An object of class \code{call}.}
#' \item{copulafam}{Prespecified Archimedean copula structure.}
#' \item{theta}{Estimated copula parameter.}
#' \item{tau.alpha}{Estimated Kendall's tau.}
#' \item{tau.theta}{Estimated Kendall's tau.}
#' \item{S_D}{Estimated survival function of the terminal event time.}
#' \item{S_k}{A list of estimated survival functions of multiple intermediate event time.}
#' \item{mar.fits}{A list of scrsurv objects.}
#' }
#' @seealso \code{\link{plot.mscr}}
#' @export
#'
#' @references Yu, Tonghui, and Xiang, Liming (2024). Association analysis of multiple intermediate events and dynamic terminal prediction.
#' \emph{Working paper.}
## #' \emph{arXiv}, \bold{1}(0): 0--0.
#' @references Yu, Tonghui, Xiang, Liming, Chen, Chixiang, and Chiou, Sy Han (2024). CopulaSCR:  .
#' \emph{Working paper.}
###' \emph{Journal of XXX}, \bold{105}(5): 1--34.
#' @examples
#' \dontrun{
#' set.seed(12345)
#' data<- simSCRmul(n=100,copulafam = "frank",K=3,
#'                  tau.alpha = .5,tau.theta=.5)
#'
#' fit<- mscr(mT1=data$T1,T2=data$T2,mevent1=data$event1,event2=data$event2,
#'            copulafam = "frank",ncore = 1)
#'
#' print(fit)
#' plot(fit$mar.fits[[1]])
#' plot(fit,type="msurv")
#' # if (!require("BiocManager", quietly = TRUE))
#' #   install.packages("BiocManager")
#' # BiocManager::install("graph")
#' # BiocManager::install("Rgraphviz")
#' plot(fit,type="data")
#' }
#'
#'
mscr<- function(mT1,T2,mevent1,event2,solver = c("optimize", "BBspg"),
                copulafam=c("frank","clayton","joe","gumbel","amh"),
                nsim=500,ncore= 4,tol=0.01,a=NULL,b = NULL,
                positive=TRUE,msurv.method= "JFKC",lower=0.02,upper=0.96){
  Call <- match.call()
  copulafam<- match.arg(copulafam)
  solver<- match.arg(solver)
  indx <- match(c("mT1","T2","mevent1","event2"), names(Call), nomatch=0)

  if (indx[1]==0) stop("a mT1 argument is required")
  if (indx[2]==0) stop("a T2 argument is required")
  if (indx[3]==0) stop("a mevent1 argument is required")
  if (indx[4]==0) stop("a event2 argument is required")

  if (!(is.matrix(mT1) | is.data.frame(mT1))) {
    stop("mT1 object must be a matrix or a data.frame")
  }

  if (!(is.matrix(mevent1) | is.data.frame(mevent1))) {
    stop("mevent1 object must be a matrix or a data.frame.")
  }

  if (!(is.matrix(T2) | is.vector(T2))) {
    stop("T2 object must be a matrix or a vector.")
  }

  if (!(is.matrix(event2) | is.vector(event2))) {
    stop("event2 object must be a matrix or a vector.")
  }

  T2<- as.vector(T2)
  event2<- as.vector(event2)

  if(ncol(mT1)!=ncol(mevent1)){
    stop("mT1 and mevent1 shall have same number of columns and rows.")
  }

  if(nrow(mT1)!=nrow(mevent1)){
    stop("mT1 and mevent1 shall have same number of columns and rows.")
  }

  if(nrow(mT1)!=length(T2)){
    stop("The length of T2 shall equal the number of rows of mT1.")
  }

  if(length(event2)!=length(T2)){
    stop("T2 and event2 shall have same length.")
  }

  K<- ncol(mT1)

  if(K==1){
    stop("The number of columns in mT1 and mevent1 shall >1.")
  }

  exact<- TRUE
  smooth<- FALSE
  jumpsize<- FALSE

  t1fits<- list()

  if(ncore>1&!(K<=3&sum(1-event2)<10)){
    Cores <- parallel::makeCluster(min(parallel::detectCores(),ncore,2^K))
    doParallel::registerDoParallel(Cores)
    parallel<- TRUE
  }else{
    parallel<- FALSE
  }


  if(isTRUE(parallel)){
    t1fits<- foreach(k = 1:ncol(mT1),
            .export = c("scrassonp","scrassonp0","assonpfit","Caltau","scrnpGEE",
                        "npCHRfit","prepformula","prepData0","prepData",
                        "copulabd_pos","copulabd"))%dopar%{
                          sdata<- data.frame(T1 = mT1[,k],T2=T2,event1=mevent1[,k],event2=event2)
                          fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
                                              t2.formula=Surv(T2, event2) ~ 1,
                                              data = sdata,copulafam=copulafam, B=0,
                                              se = FALSE,a = a,b = b,tol=1e-5,positive = positive)
                          scrsurv(fit = fitasso, method= msurv.method,
                                  surv2km=TRUE,conf.int =FALSE)
                          }
  }else{
    for(k in 1:ncol(mT1)){
      sdata<- data.frame(T1 = mT1[,k],T2=T2,event1=mevent1[,k],event2=event2)
      fitasso<- scrassonp(t1.formula=Surv(T1, event1) ~ 1,
                          t2.formula=Surv(T2, event2) ~ 1,
                          data = sdata,copulafam=copulafam,
                          se = FALSE,a = a,b = b,tol=1e-5,positive = positive)

      # cat("tau.est =", fitasso$tau,"\n")
      t1fits[[k]]<- scrsurv(fit = fitasso, method= msurv.method,
                            surv2km=TRUE,conf.int =FALSE)
      rm(fitasso,sdata)
    }
  }


  # S_D<- approxfun(x= c(0,t1fits[[1]]$t2.surv$time),y=c(1,t1fits[[1]]$t2.surv$surv),
  #                 yleft = 1,yright = min(t1fits[[1]]$t2.surv$surv),method = "linear",
  #                 f = 0,ties = mean,na.rm = TRUE)
  # S_k<- lapply(1:K,function(k)
  #   approxfun(x=c(0,t1fits[[k]]$t1.surv$time),y=c(1,t1fits[[k]]$t1.surv$surv),
  #             yleft = 1,yright = min(t1fits[[k]]$t1.surv$surv),method = "linear",
  #             f = 0,ties = mean,na.rm = TRUE))


  S_D<- stepfun(x= t1fits[[1]]$t2.surv$time,y=c(1,t1fits[[1]]$t2.surv$surv),
                right = FALSE,ties = max)
  S_k<- lapply(1:K,function(k)
    stepfun(x=t1fits[[k]]$t1.surv$time,y=c(1,t1fits[[k]]$t1.surv$surv),
            right = FALSE,ties = max))

  names(t1fits)<- names(S_k)<- colnames(mT1)

  theta<- sapply(t1fits,function(x)x$copulaparam)
  Sdev_k<- lapply(1:K,function(k)
    survnumdiff(t1fits[[k]]$t1.surv$time,t1fits[[k]]$t1.surv$surv,
                smooth=smooth,jumpsize=jumpsize,exact=exact))
  Sdev_D<- survnumdiff(t1fits[[1]]$t2.surv$time,t1fits[[1]]$t2.surv$surv,
                       smooth=smooth,jumpsize=jumpsize,exact=exact)
  # Sdev_k<- lapply(1:K,function(k)
  #   stepfun(x=t1fits[[k]]$t1.surv$time,
  #           y=c(1,diff(c(1,t1fits[[k]]$t1.surv$surv))),
  #           right = FALSE))
  # Sdev_D<- stepfun(x=t1fits[[1]]$t2.surv$time,
  #                  y=c(1,diff(c(1,t1fits[[1]]$t2.surv$surv))),
  #                  right = FALSE)


  alphafit<- mscrassofit(mT1 = mT1,T2 = T2,mevent1 = mevent1,event2 = event2,
                         S_D = S_D,S_k = S_k,Sdev_D = Sdev_D,Sdev_k = Sdev_k,
                         theta = theta,copulafam = copulafam,nsim=nsim,
                         parallel= parallel,tol=tol,lower=lower,upper=upper,solver=solver)


  if(isTRUE(parallel)){
    parallel::stopCluster(Cores)
  }

  result<- alphafit
  result$alpha.int<- c(lower,upper)
  result$call<- Call
  result$copulafam<- copulafam
  result$theta<- theta
  result$tau.alpha<- Caltau(copulafam,alphafit$alpha)
  result$tau.theta<- sapply(theta,Caltau,copulafam =copulafam)
  result$S_D<- S_D
  result$S_k<- S_k
  result$mar.fits<- t1fits
  class(result)<- "mscr"
  return(result)

}

#' @noRd
mscrassofit<- function(mT1,T2,mevent1,event2,S_D,S_k,Sdev_D,Sdev_k,
                       theta,copulafam,nsim=100,parallel= FALSE,tol=1e-1,
                       lower=0.02,upper=0.98,solver="optimize"){

  n<- length(T2)
  K<- ncol(mT1)
  di<- apply(mevent1,1,sum,na.rm = TRUE)

  # if(ncore>1&!(K<=3&sum(1-event2)<10)){
  #   Cores <- parallel::makeCluster(min(parallel::detectCores(),ncore,2^max(di)))
  #   doParallel::registerDoParallel(Cores)
  #   parallel<- TRUE
  # }else{
  #   parallel<- FALSE
  # }




  SD.Y<- sapply(T2, S_D)
  Sdev.Y<- sapply(T2, Sdev_D)

  Sk.T<- lapply(1:K, function(k)sapply(mT1[,k], S_k[[k]]))
  Sdevk.t<- lapply(1:K, function(k)sapply(mT1[,k], Sdev_k[[k]]))

  S_k_D<- lapply(1:K, function(k)sapply(1:n,function(i)
    Calcondsurv_D(surv1=Sk.T[[k]][i],surv2=SD.Y[i],
                  copulafam=copulafam,copulaparam=theta[k])))
  Sdev_k_D<- lapply(1:K, function(k)sapply(1:n,function(i)
    dev_Copula(copulafam = copulafam,param = theta[k],
               p1=Sk.T[[k]][i],p2=SD.Y[i],mode = "12"))*Sdevk.t[[k]]*
      (mT1[,k]<=T2))



  jps2<- sort(unique(c(0,T2,2*T2)),decreasing = FALSE)
  jps1<- lapply(1:K, function(k)sort(unique(c(0,mT1[,k]))))
  SD.t<- sapply(jps2, S_D)
  jps2<- jps2[c(1,which(diff(SD.t)<0)+1)]
  SD.t<- SD.t[c(1,which(diff(SD.t)<0)+1)]



  pos1<- which(event2==1)
  pos0<- which(event2==0)
  interval = copulabd_mT1(copulafam =copulafam,lower=lower,upper=upper )

  if(solver == "optimize"){
  optfit<- stats::optimize(f = function(alpha){
    -Calplik(alpha,n = n,K = K,di = di,pos0 = pos0,pos1 = pos1,jps1 = jps1,
             jps2 = jps2,mT1 = mT1,T2 = T2,mevent1 = mevent1,Sdev.Y = Sdev.Y,
             S_k_D = S_k_D,Sdev_k_D = Sdev_k_D,copulafam = copulafam,
             theta = theta,SD.t = SD.t,S_k = S_k,Sdev_k = Sdev_k,
             nsim=nsim,parallel=parallel)},tol = tol,
    interval = interval)
  alpha<- optfit$minimum
  logLik<-  -optfit$objective ### loglik/n
}
  # optfit<- optim(f = function(alpha){
  #   -Calplik(alpha,n = n,K = K,di = di,pos0 = pos0,pos1 = pos1,jps1 = jps1,
  #            jps2 = jps2,mT1 = mT1,T2 = T2,mevent1 = mevent1,Sdev.Y = Sdev.Y,
  #            S_k_D = S_k_D,Sdev_k_D = Sdev_k_D,copulafam = copulafam,
  #            theta = theta,SD.t = SD.t,S_k = S_k,Sdev_k = Sdev_k,
  #            nsim=nsim,parallel=parallel)},control = list(trace=TRUE),
  #   par = median(interval),method = "Nelder-Mead",hessian = FALSE)
  #
  # alpha<- optfit$par
  # logLik<-  -optfit$value ### loglik/n

  # as<- sapply(seq(0.01,0.95,length.out=40),Calitau,copulafam=copulafam)
  # fs<- sapply(as,f) ## convex function with one minimum
  # Caltau(copulafam,as[which.min(fs)])

  if(solver == "BBspg"){
    optfit<- BB::spg(par=0.8*interval[1]+0.2*interval[2],fn = function(alpha){
      -Calplik(alpha,n = n,K = K,di = di,pos0 = pos0,pos1 = pos1,jps1 = jps1,
               jps2 = jps2,mT1 = mT1,T2 = T2,mevent1 = mevent1,Sdev.Y = Sdev.Y,
               S_k_D = S_k_D,Sdev_k_D = Sdev_k_D,copulafam = copulafam,
               theta = theta,SD.t = SD.t,S_k = S_k,Sdev_k = Sdev_k,
               nsim=nsim,parallel=parallel)},
      lower= interval[1],upper=interval[2],quiet = FALSE,
      alertConvergence=TRUE,method = 3,
      control = list(maxit=500,ftol=tol,maxfeval=500,maximize = FALSE,eps=1e-3,
                     checkGrad=FALSE))
    alpha<- optfit$par
    logLik<-  -optfit$value ### loglik/n
  }
  # if(isTRUE(parallel)){
  #   parallel::stopCluster(Cores)
  # }
  return(list(alpha=alpha, logLik.alpha = logLik))


}



#' @noRd
copulabd_mT1<- function(copulafam,lower=0.02,upper=0.98){
  # lb<- switch(copulafam,
  #             amh = -1,
  #             clayton = -1,
  #             frank = 0.001,
  #             gumbel = 1,
  #             joe = 1)
  # ub<- switch(copulafam,
  #             amh = 1,
  #             clayton = 400,
  #             frank = 400,
  #             gumbel = 400,
  #             joe = 400)
  # return(c(lb,ub))
  return(c(Calitau(copulafam,lower),Calitau(copulafam,upper)))
}

#' @noRd
Calplik<- function(alpha,n,K,di,pos0,pos1,jps1,jps2,mT1,T2,mevent1,
                   Sdev.Y,S_k_D,Sdev_k_D,copulafam,theta,SD.t,S_k,
                   Sdev_k,nsim=500,parallel=FALSE){

  xsim<- lpxsim(n=nsim,copulafam = copulafam,param=alpha)

  L<- rep(1,n) # make 0^0=1
  L1<- Calplik1(alpha = alpha,K = K,mevent1 = mevent1,di = di,Sdev.Y = Sdev.Y,
                S_k_D = S_k_D,Sdev_k_D = Sdev_k_D,copulafam = copulafam,xsim = xsim)
  L[pos1]<- L1[pos1]

  if(length(pos0)>0){
    L[pos0]<- Calplik2(alpha = alpha,mT1 = mT1,mevent1 = mevent1,T2 = T2,
                       pos0 = pos0,di = di,jps1 = jps1,jps2 = jps2,
                       SD.t = SD.t,S_k = S_k, Sdev_k = Sdev_k,
                       copulafam=copulafam,theta = theta,
                       xsim = xsim,parallel=parallel)
  }
  # L<- L[L< max(1,IQR(L)*1.5+quantile(L,0.75))]
  # L<- L[L>.Machine$double.xmin&L<1.2]
  # L<- L[L>.Machine$double.xmin&L<=quantile(L,0.95)]
  # plot(L)
  loglik<- log(L)
  loglik[is.infinite(loglik)]<- NA
  # loglik[loglik>1]<- NA ###???
  loglik<- ifelse(loglik>10,10,loglik)
  return(ifelse(sum(!is.na(loglik))>1,mean(loglik[!is.na(loglik)]),-.Machine$double.xmax))
  # return(ifelse(sum(!is.na(loglik))>1,sum(loglik,na.rm = TRUE),-.Machine$double.xmax))
}

#' @noRd
Calplik1<- function(alpha,K,mevent1,di,Sdev.Y,S_k_D,Sdev_k_D,copulafam,xsim){

  # minval<- do.call(c,S_k_D)
  # minval<- min(minval[minval!=0])
  # minval<- min(minval,1e-15)
  phi<- sapply(1:K,function(k)
    archmCopulaLink(copulafam=copulafam,param=alpha, p=S_k_D[[k]]))
  phiallK<- apply(phi,1,sum,na.rm = TRUE)

  # psi_d<- sapply(1:length(di),function(i)
  #   psialpha_dev(u=phiallK[i],d=di[i],copulafam=copulafam,param=alpha))
  # pdev<- lapply(sort(unique(di)),function(dd){
  #   if(dd==0){
  #     function(y)archmCopulaLink_inv(copulafam=copulafam,param=alpha, y)
  #   }else{
  #     psialpha_dev(d=dd,copulafam=copulafam,param=alpha)
  #     }})
  # psi_d<- sapply(1:length(di),function(i){
  #   dd<- which(sort(unique(di))==di[i])
  #   if(di[i]==0){
  #     pdev[[dd]](phiallK[i])
  #   }else{
  #     param<- alpha
  #     y<- phiallK[i]
  #     eval(pdev[[dd]])
  #   }  })



  psi_d<- sapply(1:length(di),function(i){
    mean(exp(-xsim*phiallK[i])*(-xsim)^di[i],na.rm = TRUE)})


  # phiprod<- sapply(1:K,function(k)
  #   archmCopulaLink_dev(copulafam=copulafam,param=alpha, S_k_D[[k]])*
  # pmin(Sdev_k_D[[k]], -.Machine$double.xmin))

  phiprod<- sapply(1:K,function(k)
    archmCopulaLink_dev(copulafam=copulafam,param=alpha, S_k_D[[k]])*
      Sdev_k_D[[k]])

  # phiprod<- sapply(1:K,function(k)
  #   archmCopulaLink_dev(copulafam=copulafam,param=alpha, S_k_D[[k]]) )
  phiprod[mevent1==0]<- 1
  phiprod<- apply(phiprod,1,prod,na.rm = TRUE)

  L1<- psi_d*phiprod*Sdev.Y*(-1)^{1+di}
  # L1<- psi_d*phiprod*pmin(- .Machine$double.xmin,Sdev.Y)*(-1)^{1+di}
  # L1[psi_d==0]<- 0
  # L1[phiprod==0]<- 0
  # L1[Sdev.Y==0]<- 0
  L1[is.na(L1)]<- 0
  return(L1)

}

#' @import foreach
#' @noRd
Calplik2<- function(alpha,mT1,mevent1,T2,pos0,di,jps1,jps2,SD.t,S_k,Sdev_k,
                    copulafam,theta,xsim,parallel=FALSE){
  delta0pos<- lapply(pos0,function(xx)which(mevent1[xx,]==0))

  delta1pos<- lapply(pos0,function(xx)which(mevent1[xx,]==1))
  di2<- di[pos0]


  L2<- c()
  for(k in seq(length(pos0))){
    loc0.sets<- list(c())
    if(length(delta0pos[[k]])>1){
      loc0.sets[2:2^length(delta0pos[[k]])]<- get_subsets(delta0pos[[k]])
    }else if(length(delta0pos[[k]])==1){
      loc0.sets[2]<- delta0pos[[k]]
    }
    if(isTRUE(parallel)&length(loc0.sets)>1){
      L2i<- foreach(vec = loc0.sets,.combine = c,
                    .export = c("L2_Jint","Calcondsurv_D","archmCopulaLink",
                                "archmCopulaLink_dev","archmCopulaLink_inv",
                                "dev_Copula"))%dopar%{
                                  L2_Jint(alpha = alpha,loc0=vec,loc1=delta1pos[[k]],
                                          t1obs=mT1[pos0[k],],t2obs=T2[pos0[k]],S_k = S_k,
                                          Sdev_k = Sdev_k,jps1 = jps1,jps2=jps2,SD.t = SD.t,
                                          theta = theta,copulafam = copulafam,xsim=xsim)}
    }else{
      L2i<- sapply(loc0.sets,function(vec)
        L2_Jint(alpha = alpha,loc0=vec,loc1=delta1pos[[k]],
                t1obs=mT1[pos0[k],],t2obs=T2[pos0[k]],S_k = S_k,
                Sdev_k = Sdev_k,jps1 = jps1,jps2=jps2,SD.t = SD.t,
                theta = theta,copulafam = copulafam,xsim=xsim))
    }


    L2i<- sum(L2i,na.rm = TRUE)
    L2<- c(L2,L2i)
  }


  return(L2)

}

#' @noRd
L2_Jint<- function(alpha,loc0,loc1,t1obs,t2obs,S_k,Sdev_k,jps1,
                   jps2,SD.t,theta,copulafam,xsim){
  jps2_2<- jps2[jps2>=t2obs]
  ny<- length(jps2_2)
  if(ny==0){return(0)}
  SD.t2<- SD.t[jps2>=t2obs]
  dd<- length(loc1)

  loc00<- setdiff(1:length(t1obs),union(loc0,loc1))

  if(dd!=0&length(loc0)>=1){

    Sk.T1<- sapply(loc1, function(k)S_k[[k]](t1obs[k]))
    # Sdevk.t1<- sapply(loc1, function(k)
    #   pmin(- .Machine$double.xmin,Sdev_k[[k]](t1obs[k])))
    # Sdevk.t1<- rep(-1,length(loc1))
    Sdevk.t1<- sapply(loc1, function(k)Sdev_k[[k]](t1obs[k]))
    S_k_Dt1<- lapply(1:dd, function(k)sapply(1:ny,function(i)
      Calcondsurv_D(surv1=Sk.T1[k],surv2=SD.t2[i],
                    copulafam=copulafam,copulaparam=theta[loc1[k]])))

    phi1.t<- lapply(1:dd,function(k)
      archmCopulaLink(copulafam=copulafam,param=alpha, p=S_k_Dt1[[k]]))
    phi1.t<- Reduce("+",phi1.t)


    Sdev_k_D1<- lapply(1:dd, function(k)sapply(1:ny,function(i)
      dev_Copula(copulafam = copulafam,param = theta[loc1[k]],
                 p1=Sk.T1[k],p2=SD.t2[i],mode = "12"))*Sdevk.t1[k])
    phiprod.t<- lapply(1:dd,function(k)
      archmCopulaLink_dev(copulafam=copulafam,param=alpha, js = S_k_Dt1[[k]])*
        Sdev_k_D1[[k]] )
    phiprod.t<- Reduce("*",phiprod.t)


    Sk.t<- lapply(loc0, function(k)sapply(jps2_2, S_k[[k]]))
    Sk.Y<- sapply(loc0, function(k)S_k[[k]](t2obs))

    S_k_D0.tt<- lapply(1:length(loc0), function(k)sapply(1:ny,function(i)
      Calcondsurv_D(surv1=Sk.t[[k]][i],surv2=SD.t2[i],
                    copulafam=copulafam,copulaparam=theta[loc0[k]])))
    S_k_D0.Yt<- lapply(1:length(loc0), function(k)
      Calcondsurv_D(surv1=Sk.Y[k],surv2=SD.t2,
                    copulafam=copulafam,copulaparam=theta[loc0[k]]))


    phi0.tt<- lapply(1:length(loc0),function(k)
      sapply(S_k_D0.tt[[k]],function(xx)
        archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))

    phi0.Yt<- lapply(1:length(loc0),function(k)
      sapply(S_k_D0.Yt[[k]],function(xx)
        archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))

    if(length(loc00)!=0){
      Sk.t00<- lapply(loc00, function(k)sapply(jps2_2, S_k[[k]]))
      S_k_D00.tt<- lapply(1:length(loc00), function(k)sapply(1:ny,function(i)
        Calcondsurv_D(surv1=Sk.t00[[k]][i],surv2=SD.t2[i],
                      copulafam=copulafam,copulaparam=theta[loc00[k]])))
      phi00.tt<- lapply(1:length(loc00),function(k)
        sapply(S_k_D00.tt[[k]],function(xx)
          archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))
      phi00.tt<- Reduce("+",phi00.tt)
    }else{
      phi00.tt<- rep(0,ny)
    }



    phi0.x.t<- lapply(1:ny, function(ti)
      Reduce("*",lapply(1:length(loc0),function(k)
        exp(-phi0.tt[[k]][ti]*xsim)-exp(-phi0.Yt[[k]][ti]*xsim))))

    phi.t<- sapply(1:ny, function(ti)
      mean(exp(-xsim*phi1.t[ti])*exp(-xsim*phi00.tt[ti])* phi0.x.t[[ti]]*(-xsim)^dd,na.rm = TRUE))


    Js<- sum(diff(c(SD.t2,0))*phi.t* phiprod.t,na.rm = TRUE)*(-1)^{dd+1+length(loc0)}



  }else if(dd!=0&length(loc0)==0){
    Sk.T1<- sapply(loc1, function(k)S_k[[k]](t1obs[k]))
    Sdevk.t1<- sapply(loc1, function(k)Sdev_k[[k]](t1obs[k]))
    # Sdevk.t1<- sapply(loc1, function(k)
    # pmin(- .Machine$double.xmin,Sdev_k[[k]](t1obs[k])))
    # Sdevk.t1<- rep(-1,length(loc1))
    S_k_Dt1<- lapply(1:dd, function(k)sapply(1:ny,function(i)
      Calcondsurv_D(surv1=Sk.T1[k],surv2=SD.t2[i],
                    copulafam=copulafam,copulaparam=theta[loc1[k]])))

    phi1.t<- lapply(1:dd,function(k)
      archmCopulaLink(copulafam=copulafam,param=alpha, p=S_k_Dt1[[k]]))
    phi1.t<- Reduce("+",phi1.t)

    Sdev_k_D1<- lapply(1:dd, function(k)sapply(1:ny,function(i)
      dev_Copula(copulafam = copulafam,param = theta[loc1[k]],
                 p1=Sk.T1[k],p2=SD.t2[i],mode = "12"))*Sdevk.t1[k])
    phiprod.t<- lapply(1:dd,function(k)
      archmCopulaLink_dev(copulafam=copulafam,param=alpha, js = S_k_Dt1[[k]])*
        Sdev_k_D1[[k]] )
    phiprod.t<- Reduce("*",phiprod.t)

    if(length(loc00)!=0){
      Sk.t00<- lapply(loc00, function(k)sapply(jps2_2, S_k[[k]]))
      S_k_D00.tt<- lapply(1:length(loc00), function(k)sapply(1:ny,function(i)
        Calcondsurv_D(surv1=Sk.t00[[k]][i],surv2=SD.t2[i],
                      copulafam=copulafam,copulaparam=theta[loc00[k]])))
      phi00.tt<- lapply(1:length(loc00),function(k)
        sapply(S_k_D00.tt[[k]],function(xx)
          archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))
      phi00.tt<- Reduce("+",phi00.tt)
    }else{
      phi00.tt<- rep(0,ny)
    }

    # pd<- psialpha_dev(d=dd,copulafam=copulafam,param=alpha)
    #
    # y<- phi1.t+phi00.tt
    # param<- alpha
    # psi_d<- eval(pd)
    #
    # Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t,na.rm = TRUE)*(-1)^{dd+1}



    psi_d<- sapply(1:ny, function(ti)
      mean(exp(-xsim*phi1.t[ti])*exp(-xsim*phi00.tt[ti])*(-xsim)^dd,na.rm = TRUE))

    Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t,na.rm = TRUE)*(-1)^{dd+1}
    if(!is.na(Js)){
      if(Js>1){
        psi_d<- sapply(1:ny, function(ti)
          median(exp(-xsim*phi1.t[ti])*exp(-xsim*phi00.tt[ti])*(-xsim)^dd,na.rm = TRUE))

        Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t,na.rm = TRUE)*(-1)^{dd+1}
      }
    }

  }else if(dd==0&length(loc0)==0){
    if(length(loc00)!=0){
      Sk.t00<- lapply(loc00, function(k)sapply(jps2_2, S_k[[k]]))
      S_k_D00.tt<- lapply(1:length(loc00), function(k)sapply(1:ny,function(i)
        Calcondsurv_D(surv1=Sk.t00[[k]][i],surv2=SD.t2[i],
                      copulafam=copulafam,copulaparam=theta[loc00[k]])))
      phi00.tt<- lapply(1:length(loc00),function(k)
        sapply(S_k_D00.tt[[k]],function(xx)
          archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))
      phi00.tt<- Reduce("+",phi00.tt)
    }else{
      phi00.tt<- rep(0,ny)
    }


    psi_d<- archmCopulaLink_inv(copulafam=copulafam,param=alpha, y = phi00.tt)


    Js<- sum(diff(c(SD.t2,0))*psi_d,na.rm = TRUE)*(-1)



  }else if(dd==0&length(loc0)>=1){

    Sk.t<- lapply(loc0, function(k)sapply(jps2_2, S_k[[k]]))
    Sk.Y<- sapply(loc0, function(k)S_k[[k]](t2obs))

    S_k_D0.tt<- lapply(1:length(loc0), function(k)sapply(1:ny,function(i)
      Calcondsurv_D(surv1=Sk.t[[k]][i],surv2=SD.t2[i],
                    copulafam=copulafam,copulaparam=theta[loc0[k]])))
    S_k_D0.Yt<- lapply(1:length(loc0), function(k)
      Calcondsurv_D(surv1=Sk.Y[k],surv2=SD.t2,
                    copulafam=copulafam,copulaparam=theta[loc0[k]]))


    phi0.tt<- lapply(1:length(loc0),function(k)
      sapply(S_k_D0.tt[[k]],function(xx)
        archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))

    phi0.Yt<- lapply(1:length(loc0),function(k)
      sapply(S_k_D0.Yt[[k]],function(xx)
        archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))

    if(length(loc00)!=0){
      Sk.t00<- lapply(loc00, function(k)sapply(jps2_2, S_k[[k]]))
      S_k_D00.tt<- lapply(1:length(loc00), function(k)sapply(1:ny,function(i)
        Calcondsurv_D(surv1=Sk.t00[[k]][i],surv2=SD.t2[i],
                      copulafam=copulafam,copulaparam=theta[loc00[k]])))
      phi00.tt<- lapply(1:length(loc00),function(k)
        sapply(S_k_D00.tt[[k]],function(xx)
          archmCopulaLink(copulafam=copulafam,param=alpha, p=xx)))
      phi00.tt<- Reduce("+",phi00.tt)
    }else{
      phi00.tt<- rep(0,ny)
    }


    phi0.x.t<- lapply(1:ny, function(ti)
      Reduce("*",lapply(1:length(loc0),function(k)
        exp(-phi0.tt[[k]][ti]*xsim)-exp(-phi0.Yt[[k]][ti]*xsim))))

    phi.t<- sapply(1:ny, function(ti)
      mean(exp(-xsim*phi00.tt[ti])*phi0.x.t[[ti]],na.rm = TRUE))


    Js<- sum(diff(c(SD.t2,0))*phi.t,na.rm = TRUE)*(-1)^{1+length(loc0)}



  }

  return(Js)
}

#' @noRd
lpxsim<- function(n,copulafam,param){
  if(copulafam=="amh"){
    xsim<- copula::copAMH@V0(n = n, theta = param) # geometric
  }
  if(copulafam=="clayton"){
    xsim<-copula::copClayton@V0(n = n, theta = param) # gamma
  }
  if(copulafam=="frank"){
    xsim<-copula::copFrank@V0(n = n, theta = param) #logarithmic series
  }
  if(copulafam=="gumbel"){
    xsim<-copula::rstable1(n = n,alpha = 1/param,beta = 1,
                           gamma = (cos(pi/2/param))^{param},delta = 0,pm = 1)
    # positive stable
  }
  if(copulafam=="joe"){
    xsim<-copula::copJoe@V0(n = n, theta = param) # Sibuya latent (power series)
  }

  return(xsim)
}



# psialpha_dev<- function(d,copulafam,param){
#   # psialpha_dev<- function(u,d,copulafam,param){
#   # if(copulafam=="amh"){
#   #   phi_inv<-  "(param-1)/(param - exp(y))"
#   # }
#   # if(copulafam=="clayton"){
#   #   phi_inv<- "(param*y+1)^{-1/param}"
#   # }
#   # if(copulafam=="frank"){
#   #   phi_inv<- "-1/param*log(1+(exp(-param)-1)*exp(y))"
#   # }
#   # if(copulafam=="gumbel"){
#   #   phi_inv<- "exp(-y^(1/param))"
#   # }
#   # if(copulafam=="joe"){
#   #   phi_inv<-"1- (1-exp(-y))^(1/param)"
#   # }
#   # calculus::derivative(f = phi_inv, var = c(y = u),
#   #                      order = d,params = list(param=param))
#
#
#   if(copulafam=="amh"){
#     phi_inv<- expression((param-1)/(param - exp(y)))
#   }
#   if(copulafam=="clayton"){
#     phi_inv<- expression((param*y+1)^{-1/param})
#   }
#   if(copulafam=="frank"){
#     phi_inv<- expression(-1/param*log(1+(exp(-param)-1)*exp(y)))
#   }
#   if(copulafam=="gumbel"){
#     phi_inv<- expression(exp(-y^(1/param)))
#   }
#   if(copulafam=="joe"){
#     phi_inv<-expression(1- (1-exp(-y))^(1/param))
#   }
#
#   df<-DD(phi_inv,"y",d)
#
#
#   df
# }
# DD <- function(expr, name, order = 1,...) {
#   if(order < 1) stop("'order' must be >= 1")
#   if(order == 1) D(expr, name,...)
#   else DD(D(expr, name,...), name, order - 1,...)
# }
#' @noRd
get_subsets<- function(vec){
  vec<- as.numeric(unique(vec))
  unlist(lapply(1:length(vec),
                combn,
                x = vec,
                simplify = FALSE),
         recursive = FALSE)

}
#' @noRd
Calcondsurv_D<- function(surv1,surv2,copulafam,copulaparam){
  if(length(surv1)>1&length(surv2)>1){
    out<- lapply(surv2, function(xx)
      dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =xx,mode = "2"))
  }else{
    out<- dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =surv2,mode = "2")
  }
  return(out)

}

#' @noRd
survnumdiff<- function(x,y,smooth,jumpsize,span = 0.75,exact=FALSE){

  if(!isTRUE(exact)){

    # da<- data.frame(x=x.new,y=y.new)
    if(isTRUE(jumpsize)){
      y<- c(1,y)
      x<- c(0,x)
      # x.new<- x[!is.na(y)]
      # y.new<- y[!is.na(y)]
      ypos<- unique(y[!is.na(y)])
      ypos<- sapply(ypos,function(yy)which(y==yy)[1])
      x.new<- x[ypos]
      y.new<- y[ypos]
      sdev<-  diff(c(1,y.new),lag = 1)
    }else{
      y<- c(1,y)
      x<- c(0,x)
      ypos<- unique(y[!is.na(y)])
      ypos<- sapply(ypos,function(yy)which(y==yy)[1])
      x.new<- x[ypos]
      y.new<- y[ypos]
      sdev<-  diff(y.new,lag = 2)/diff(x.new,lag = 2)
      sdev<- c((y.new[2]-y.new[1])/(x.new[2]-x.new[1]),sdev)
      sdev<- c(sdev,(y.new[length(y.new)]-y.new[length(y.new)-1])/
                 (x.new[length(x.new)]-x.new[length(x.new)-1]))
      x.new<- x.new[!is.na(sdev)&!is.infinite(sdev)]
      sdev<- sdev[!is.na(sdev)&!is.infinite(sdev)]

    }
    if(isTRUE(smooth)){
      sdev<- predict(stats::loess(sdev~x.new,span=span),x.new)
    }

    sdev<- stats::approxfun(x=x.new,y=sdev,method = "linear",
                     yleft = 0,yright = 0,f=0,na.rm = TRUE,ties = mean)
  }

  if(isTRUE(exact)){
    y<- c(1,y)
    x<- c(0,x)
    x.new<- x[!is.na(y)]
    y.new<- y[!is.na(y)]

    xpos<- unique(x.new)
    xpos<- sapply(xpos,function(xx)which(x.new==xx)[1])
    x.new<- x.new[xpos]
    y.new<- y.new[xpos]

    # sdev0<-  diff(c(1,y.new),lag = 1)


    # sdev<- function(xx){
    #   xpos<- match(xx,x.new,nomatch = FALSE)
    #   yy<- numeric(length = length(xx))
    #   yy[xpos!=0]<- sdev0[xpos[xpos!=0]]
    #   yy
    # }
    # sdev<- approxfun(x=x.new,y=sdev0,method = "constant",
    #                  yleft = 0,yright = 0,f=0,na.rm = TRUE,ties = mean)


    sdev<- stats::stepfun(x=x.new,y=c(0,diff(c(1,y.new))), right = FALSE,ties = mean )
  }



  return(sdev)

}
