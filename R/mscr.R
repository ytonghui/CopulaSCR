# library(calculus)
library(doParallel)
source("CalCopula.R")
library(survival)
##### main functions #############
mscrdyn<- function(mT1,T2,mevent1,event2,
                   copulafam=c("frank","clayton","joe","gumbel","amh"),
                   nsim=500,ncore= 1,tol=0.01,a=NULL,b = NULL,
                   smooth= FALSE,jumpsize=TRUE,exact= TRUE,
                   positive=TRUE,msurv.method= "JFKC",lower=0.02,upper=0.98){
  Call <- match.call()
  copulafam<- match.arg(copulafam)
  
  indx <- match(c("mT1","T2","mevent1","event2"), names(Call), nomatch=0)
  
  if (indx[1]==0) stop("a mT1 argument is required")
  if (indx[2]==0) stop("a T2 argument is required")
  if (indx[3]==0) stop("a mevent1 argument is required")
  if (indx[4]==0) stop("a event2 argument is required")
  
  if (!(is.matrix(mT1) | is.data.frame(mT1))) {
    stop("mT1 object must be a matrix or a data.frame")
  }
  
  if (!(is.matrix(mevent1) | is.data.frame(mevent1))) {
    stop("mevent1 object must be a matrix or a data.frame")
  }
  
 
  
  if(ncol(mT1)!=ncol(mevent1)){
    stop("mT1 and mevent1 shall have same number of columns and rows.")
  }
  
  K<- ncol(mT1)
  
  if(K==1){
    stop("The number of columns in mT1 and mevent1 shall >1.")
  }
  
  t1fits<- list()
  for(k in 1:ncol(mT1)){
    sdata<- data.frame(T1 = mT1[,k],T2=T2,event1=mevent1[,k],event2=event2)
    fitasso<- scrassonp(sdata,t1.formula=Surv(T1, event1) ~ 1, 
                        t2.formula=Surv(T2, event2) ~ 1,copulafam=copulafam,
                        model=TRUE,se = FALSE,a = a,b = b,tol=1e-5,
                        positive = positive)
    # cat("tau.est =", fitasso$tau,"\n")
    t1fits[[k]]<- scrsurv(fit = fitasso, method= msurv.method,surv2km=TRUE,confint=FALSE)
    rm(fitasso)
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
                         ncore= ncore,tol=tol,lower=lower,upper=upper)
  result<- alphafit
  result$call<- Call
  result$copulafam<- copulafam
  result$theta<- theta
  result$tau.alpha<- Caltau(copulafam,alphafit$alpha)
  result$tau.theta<- sapply(theta,Caltau,copulafam =copulafam)
  result$S_D<- S_D
  result$S_k<- S_k
  result$t1fits<- t1fits
  class(result)<- "mscr"
  return(result)
  
}


predict.mscr<- function(fit, t1obs,times,type=c("surv", "rrms", "rmst", "qrl","qst"),
                        tu,tau=0.5,nsim=1000,smooth= FALSE,jumpsize=FALSE,exact=FALSE){
  #times	：Vector of times at which to return the estimated survival rates.
  Call <- match.call()
  indx <- match(c('fit', 't1obs','times','tu'), names(Call), nomatch=0)
  type<- match.arg(type)
  
  if (indx[1]==0) stop("a fit argument of class 'mscr' is required")
  if (indx[2]==0) stop("a t1obs argument is required")
  
  
  if(class(fit)!= "mscr"){
    stop("fit shall be 'mscr' class.")
  }
  
  if (indx[3]==0&type=="surv") {
    times<- knots(fit$S_D)
  }
  if(type!="surv"){
    times<- knots(fit$S_D)
  }
  
  
  if (indx[4]==0){
    tu<- max(knots(fit$S_D))
  }
  
  times<- times[times<= tu]
  
  if (!(is.vector(t1obs)|is.matrix(t1obs) | is.data.frame(t1obs))) {
    stop("t1obs object must be a vector, matrix or a data.frame")
  }
  
  if(is.vector(t1obs)){
    if(is.null(names(t1obs))){
      stop("Please give each element of the vector t1obs a name to specify 
         which intermediate event times are observed.")
    }else{
      if(any(!names(t1obs)%in%names(fit$S_k))){
        stop("names of the vector t1obs shall be part of colnames of 'mT1'.")
      }
      
    }
    t1vec<- TRUE
  }
  if(is.matrix(t1obs) | is.data.frame(t1obs)){
    if(is.null(colnames(t1obs))){
      stop("Please give each element of the vector t1obs a name to specify 
         which intermediate event times are observed.")
    }else{
      if(any(!colnames(t1obs)%in%names(fit$S_k))){
        stop("names of the vector t1obs shall be part of colnames of 'mT1'.")
      }
      
    }
    t1vec<- FALSE
  }
  
  
  if(t1vec){
    out<- dypredfit(fit = fit,t1obs = t1obs,type = type,times = times,
                    nsim = nsim,tu = tu,tau = tau,
                    smooth=smooth,jumpsize=jumpsize,exact=exact)
  }
  
  if(!t1vec){
    # t1obs<- t1obs[apply(is.na(t1obs),1,sum)<ncol(t1obs),]
    out<- lapply(1:nrow(t1obs),function(i)
      dypredfit(fit = fit,t1obs = t1obs[i,], 
                type = type,times = times,nsim = nsim,
                tu = tu,tau = tau,smooth=smooth,jumpsize=jumpsize,exact=exact))
    out<- list(out,t1=t1obs)
    names(out)[1]<- type
  }
  return(out)
  
  
}


predbsmscr<- function(fit, t1obs,times,type=c("surv","rrms", "rmst", "qrl","qst"),
                      tu,tau=0.5,cause=0,maxobs=TRUE){
  #times	：Vector of times at which to return the estimated survival rates.
  Call <- match.call()
  indx <- match(c('fit', 't1obs','times','tu'), names(Call), nomatch=0)
  type<- match.arg(type)
  
  if (indx[1]==0) stop("a fit argument of class 'mscr' is required")
  if (indx[2]==0) stop("a t1obs argument is required")
  
  
  if(class(fit)!= "mscr"){
    stop("fit shall be 'mscr' class.")
  }
  
  if (indx[3]==0&type=="surv") {
    times<- knots(fit$S_D)
  }
  if(type!="surv"){
    times<- knots(fit$S_D)
  }
  
  
  if (indx[4]==0){
    tu<- max(knots(fit$S_D))
  }
  
  times<- times[times<= tu]
  if (!(is.vector(t1obs)|is.matrix(t1obs) | is.data.frame(t1obs))) {
    stop("t1obs object must be a vector, matrix or a data.frame")
  }
  
  if(is.vector(t1obs)){
    if(is.null(names(t1obs))){
      stop("Please give each element of the vector t1obs a name to specify 
         which intermediate event times are observed.")
    }else{
      if(any(!names(t1obs)%in%names(fit$S_k))){
        stop("names of the vector t1obs shall be part of colnames of 'mT1'.")
      }
      
    }
    t1vec<- TRUE
  }
  if(is.matrix(t1obs) | is.data.frame(t1obs)){
    if(is.null(colnames(t1obs))){
      stop("Please give each element of the vector t1obs a name to specify 
         which intermediate event times are observed.")
    }else{
      if(any(!colnames(t1obs)%in%names(fit$S_k))){
        stop("names of the vector t1obs shall be part of colnames of 'mT1'.")
      }
      
    }
    t1vec<- FALSE
  }
  
  
  if(t1vec){
    out<- bspredfit(fit = fit,t1obs = t1obs,type = type,times = times,
                    tu = tu,tau = tau,cause = cause,maxobs=maxobs)
  }
  
  if(!t1vec){
    # t1obs<- t1obs[apply(is.na(t1obs),1,sum)<ncol(t1obs),]
    out<- lapply(1:nrow(t1obs),function(i)
      bspredfit(fit = fit,t1obs = t1obs[i,],type = type,times = times,
                tu = tu,tau = tau,cause = cause,maxobs=maxobs))
    out<- list(out,t1=t1obs)
    names(out)[1]<- type
  }
  return(out)
  
  
}

plotdyBS<- function(surv,S_D,T2,event2,times=knots(S_D),reference=TRUE,tu=NULL){
  
  times<- unique(c(0,times))
  
  if(!is.null(tu)){times<- times[times<= tu]}
  if(isTRUE(reference)){
    BS00<- plotbsBS(S_D=S_D,T2=T2,event2=event2,times=times)$BS$BS
  }
  
  
  cens.fit = prodlim::prodlim(Surv(T2,event2)~1,
                              data = data.frame(T2 = T2,event2 =event2),reverse=TRUE)
  ipcw = predict(cens.fit,times=times,level.chaos=1,mode="matrix",type="surv")
  ipcw.obs = prodlim::predictSurvIndividual(cens.fit,lag=1)
  ipcw[ipcw==0]=min(ipcw[ipcw!=0],2e-10)/2
  ipcw.obs[ipcw.obs==0]=min(ipcw.obs[ipcw.obs!=0],2e-10)/2
  
  SurvMatrix<- lapply(surv$surv,function(x)
    if(sum(!is.na(x$surv))>0){
      approx(x=x$times,y=x$surv,xout = times,method = "constant",ties = mean,
             yleft = NA,yright = min(x$surv,na.rm = TRUE),f = 0,na.rm = TRUE)$y
    }else{
      rep(NA,length(times))
    })
  
  
  SurvMatrix2<- do.call(rbind,SurvMatrix)
  BS=sapply(times, function(t){
    mean(as.numeric(T2<=t)*event2*(SurvMatrix2[,times==t]^2)/ipcw.obs,na.rm =TRUE)+
      mean(as.numeric(T2>t)*(1-SurvMatrix2[,times==t])^2/ipcw[times==t],na.rm =TRUE)
  })
  
  
  IBS=1/diff(range(times))*sum(BS*diff(c(times,max(times))),na.rm =TRUE)
  
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




###### auxillary functions ###########


dypredfit<- function(fit,t1obs,type,times,nsim,tu,tau,smooth,jumpsize,exact=FALSE){
  ## type = "surv": survival function
  ## type = "rrms": restricted residual mean survival
  ## type = "rmst": restricted residual mean survival+ landmark time rrms+max(t1obs)
  ## type = "qrl": quantile residual lifetime
  ## type = "qst": quantile survival time (qrl+ landmark time)
  
  jps<- knots(fit$S_D)
  njps<- length(jps)
  
  t1name<- names(t1obs)
  t1name<- t1name[!is.na(t1obs)]
  t1obs<- t1obs[!is.na(t1obs)]
  names(t1obs)<- t1name
  # if(length(t1obs)==0){
  #   stop("Please provide valid values for t1obs argument.")
  # }
  t1pos<- match(t1name,names(fit$S_k))
  lpos<- length(t1pos)
  
  
  if(lpos==0){
    cs<- fit$S_D(times)
    t1obs<- 0
    if(type=="surv"){
      return(list(times=times,surv = cs,t1=NA))
    }
  }
  if(lpos==1){
    cs<- predcsurv(fit=fit$t1fits[[t1pos]],t1=t1obs, t2=times, s2=t1obs, 
                   t1equal =TRUE,method="sp")$condsurv[[1]][,1]
    if(type=="surv"){
      return(list(times=times,surv = cs,t1=t1obs))
    }
  }
  if(lpos>1){
    
    dSD<- diff(fit$S_D(jps))
    Qstep<- dypred.Q(fit = fit,t1obs = t1obs,t1pos = t1pos,nsim=nsim,
                     smooth=smooth,jumpsize=jumpsize,exact=exact)
    denominator<- sum(Qstep[seq(njps-1)]*dSD*(jps[seq(njps-1)]>=max(t1obs)),na.rm = TRUE)
    
    numerator<- sapply(times,function(t)
      sum(Qstep[seq(njps-1)]*dSD*(jps[seq(njps-1)]>=t),na.rm = TRUE))
    
    cs<- rep(NA,length(times))
    cs[times>=max(t1obs)]<- numerator[times>=max(t1obs)]/denominator
    if(type=="surv"){
      return(list(times=times,surv = cs,t1=t1obs))
    }
    
    
    
  }
  
  if(type%in%c("rrms", "rmst")){
    tjp<- c(max(t1obs),times[times>=max(t1obs)])
    cs<- c(1,cs[times>=max(t1obs)])
    rrms<- sum(cs[seq(length(cs)-1)]*diff(tjp)* (tjp[seq(length(cs)-1)]<=tu),na.rm = TRUE)
  }
  if(type=="rrms"){
    return(list(rrms=rrms,tu=tu,t1=t1obs))
  }
  
  if(type=="rmst"){
    rmst<- pmin(tu,rrms+max(t1obs))
    return(list(rmst=rmst,tu=tu,t1=t1obs))
  }
  
  if(type%in%c("qrl", "qst")){
    tjp<- c(max(t1obs),times[times>=max(t1obs)])
    cs<- c(1,cs[times>=max(t1obs)])
    tau<- 1-tau
    tau<- pmin(max(cs),pmax(tau,min(cs)))
    if(length(cs)<=2){
      qrl<- NA
    }else if(all(is.na(cs[cs!=1]))){
      qrl<- NA
    }else{
      qrl<- approx(x = cs,y=tjp,xout = tau,method = "linear",ties = min)$y
    }
    
  }
  if(type=="qrl"){
    return(list(qrl=qrl-max(t1obs),tau =1-tau,t1=t1obs))
  }
  if(type=="qst"){
    return(list(qst=qrl,tau =1-tau,t1=t1obs))
  }
}

bspredfit<- function(fit,t1obs,type,times,tu,tau,cause=0,maxobs=TRUE){
  jps<- knots(fit$S_D)
  njps<- length(jps)
  if(all(is.na(t1obs))){
    tmax<- NA
  }else{
    tmax<- max(t1obs,na.rm = TRUE)
  }
  
  
  t1name<- names(t1obs)
  t1name<- t1name[!is.na(t1obs)]
  t1obs<- t1obs[!is.na(t1obs)]
  names(t1obs)<- t1name
 
  
  if(cause!=0){
    t1pos<- match(t1name,names(fit$S_k))
    t1pos<- intersect(t1pos,cause)
  }else{
      t1pos<- c()
    }
  lpos<- length(t1pos)
  
  t1obs2<- t1obs[names(fit$S_k)[t1pos]]
  
  if(lpos==0){
    cs<- fit$S_D(times)
    
    if(!is.na(tmax)|length(tmax)==1){
      cs<- cs/fit$S_D(tmax)
      cs[times<tmax]<- NA
    }else{
      tmax<- NA
    }
    
    if(type=="surv"){
      return(list(times=times,surv = cs,t1=tmax))
    }
  }
  if(lpos==1){
    if(isTRUE(maxobs)){
      cs<- predcsurv(fit=fit$t1fits[[t1pos]],t1=t1obs2, t2=times, s2=tmax, 
                     t1equal =TRUE,method="sp")$condsurv[[1]][,1]
    }else{
      cs<- predcsurv(fit=fit$t1fits[[t1pos]],t1=t1obs2, t2=times, s2=t1obs2, 
                     t1equal =TRUE,method="sp")$condsurv[[1]][,1]
    }
    
    if(type=="surv"){
      return(list(times=times,surv = cs,t1=t1obs))
    }
  }
  
  
  if(type%in%c("rrms", "rmst")){
    if(is.na(tmax)){tmax<-0;t1obs2<- 0}
    if(isTRUE(maxobs)|lpos==0){
      tjp<- c(tmax,times[times>=tmax])
      cs<- c(1,cs[times>=tmax])
    }else{
      tjp<- c(t1obs2,times[times>=t1obs2])
      cs<- c(1,cs[times>=t1obs2])
    }
    
    rrms<- sum(cs[seq(length(cs)-1)]*diff(tjp)* (tjp[seq(length(cs)-1)]<=tu),na.rm = TRUE)
    
  }
  if(type=="rrms"){
    return(list(rrms=as.numeric(rrms),tu=tu,t1=t1obs))
  }
  
  if(type=="rmst"){
    if(isTRUE(maxobs)|lpos==0){
      rmst<- rrms+tmax
    }else{
      rmst<- rrms+t1obs2
    }
    return(list(rmst=as.numeric(rmst),tu=tu,t1=t1obs))
  }
  
  if(type%in%c("qrl", "qst")){
    if(is.na(tmax)){tmax<-0;t1obs2<- 0}
    if(isTRUE(maxobs)|lpos==0){
      tjp<- c(tmax,times[times>=tmax])
      cs<- c(1,cs[times>=tmax])
    }else{
      tjp<- c(t1obs2,times[times>=t1obs2])
      cs<- c(1,cs[times>=t1obs2])
    }
    
    tau<- 1-tau
    tau<- pmin(max(cs),pmax(tau,min(cs)))
    
    if(length(cs)<=2){
      qst<- NA
    }else if(all(is.na(cs[cs!=1]))){
      qst<- NA
    }else{
      qst<- approx(x = cs,y=tjp,xout = tau,method = "linear",ties = min)$y
    }
    
  }
  if(type=="qst"){
    return(list(qst=as.numeric(qst),tau =1-tau,t1=t1obs))
    
  }
  if(type=="qrl"){
    if(isTRUE(maxobs)|lpos==0){
      qrl=qst-tmax
    }else{
      qrl=qst-t1obs2
    }
    return(list(qrl=as.numeric(qrl),tau =1-tau,t1=t1obs))
  }
}


plotbsBS<- function(S_D,T2,event2,times=knots(fit$S_D)){
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

dypred.Q<- function(fit,t1obs,t1pos,nsim=500,smooth,jumpsize,exact=FALSE){
  jps<- knots(fit$S_D)
  SD.jp<- sapply(jps, fit$S_D)
  lpos<- length(t1obs)
  njp<- length(jps)
  
  Sk.T<- sapply(1:lpos, function(k)fit$S_k[[t1pos[k]]](t1obs[k]))
  
  # if(!isTRUE(smooth)){
    # Sdev_k<- lapply(1:lpos,function(k)
    #   stepfun(x=fit$t1fits[[t1pos[k]]]$t1.surv$time,
    #           y=c(0,diff(c(1,fit$t1fits[[t1pos[k]]]$t1.surv$surv))),
    #           right = FALSE))
    
  # }else{
    # Sdev_k<- lapply(1:lpos,function(k)
    #   data.frame(x=fit$t1fits[[t1pos[k]]]$t1.surv$time,
    #              y=diff(c(1,fit$t1fits[[t1pos[k]]]$t1.surv$surv))))
    # Sdev_k<- lapply(1:lpos,function(k)Sdev_k[[k]][Sdev_k[[k]]$y!=0,])
    # Sdev_k<- lapply(1:lpos,function(k)
    #   stepfun(x=Sdev_k[[k]]$x,y=c(0,Sdev_k[[k]]$y), right = FALSE)) 
  # }
  
  Sdev_k<- lapply(1:lpos,function(k)
    survnumdiff(fit$t1fits[[t1pos[k]]]$t1.surv$time,
                fit$t1fits[[t1pos[k]]]$t1.surv$surv,
                smooth=smooth,jumpsize = jumpsize,exact=exact))
  
  Sdevk.t<- sapply(1:lpos, function(k)Sdev_k[[k]](t1obs[k]))
  
  S_k_D<- lapply(1:lpos, function(k)sapply(1:njp,function(i)
    Calcondsurv_D(surv1=Sk.T[k],surv2=SD.jp[i],
                 copulafam=fit$copulafam,copulaparam=fit$theta[t1pos[k]])))
  Sdev_k_D<- lapply(1:lpos, function(k)sapply(1:njp,function(i)
    dev_Copula(copulafam = fit$copulafam,param = fit$theta[t1pos[k]],
               p1=Sk.T[k],p2=SD.jp[i],mode = "12"))*Sdevk.t[k])
  
  xsim<- lpxsim(n=nsim,copulafam = fit$copulafam,param=fit$alpha)
  
  phi<- sapply(1:lpos,function(k)
    archmCopulaLink(copulafam=fit$copulafam,param=fit$alpha, p=S_k_D[[k]]))
  phiallK<- apply(phi,1,sum)
  
  psi_d<- sapply(1:njp,function(i){
    mean(exp(-xsim*phiallK[i])*(-xsim)^lpos,na.rm = TRUE)})
  
  
  phiprod<- sapply(1:lpos,function(k)
    archmCopulaLink_dev(copulafam=fit$copulafam,param=fit$alpha, S_k_D[[k]])*
      Sdev_k_D[[k]] ) 
  
  phiprod<- as.numeric(apply(phiprod,1,prod,na.rm = TRUE))
  
  
  Q<- psi_d*phiprod*(-1)^{lpos}
  return(Q)
}





mscrassofit<- function(mT1,T2,mevent1,event2,S_D,S_k,Sdev_D,Sdev_k,
                          theta,copulafam,nsim=100,ncore= 4,tol=1e-1,
                       lower=0.02,upper=0.98){
  
  n<- length(T2)
  K<- ncol(mT1)
  di<- apply(mevent1,1,sum)
  
  if(ncore>1&!(K<=3&sum(1-event2)<10)){
    Cores <- makeCluster(min(detectCores(),ncore,2^max(di))) 
    registerDoParallel(Cores)
    parallel<- TRUE
  }else{
    parallel<- FALSE
  }
  
  
  
  
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
  
  optfit<- optimize(f = function(alpha){
    -Calplik(alpha,n = n,K = K,di = di,pos0 = pos0,pos1 = pos1,jps1 = jps1,
            jps2 = jps2,mT1 = mT1,T2 = T2,mevent1 = mevent1,Sdev.Y = Sdev.Y,
            S_k_D = S_k_D,Sdev_k_D = Sdev_k_D,copulafam = copulafam,
            theta = theta,SD.t = SD.t,S_k = S_k,Sdev_k = Sdev_k,
            nsim=nsim,parallel=parallel)},tol = tol,
          interval = interval)
  alpha<- optfit$minimum
  logLik<-  -optfit$objective ### loglik/n
  
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
  
  
  if(isTRUE(parallel)){
    stopCluster(Cores)
  }
  return(list(alpha=alpha, logLik.alpha = logLik))
  
  
}

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
  L<- L[L>.Machine$double.xmin&L<1.2]
  # plot(L)
  loglik<- log(L)
  loglik[is.infinite(loglik)]<- NA
  # loglik[loglik>1]<- NA ###???
  
  return(ifelse(sum(!is.na(loglik))>1,mean(loglik[!is.na(loglik)]),-.Machine$double.xmax))
  # return(ifelse(sum(!is.na(loglik))>1,sum(loglik,na.rm = TRUE),-.Machine$double.xmax))
}

Calplik1<- function(alpha,K,mevent1,di,Sdev.Y,S_k_D,Sdev_k_D,copulafam,xsim){
  
  # minval<- do.call(c,S_k_D)
  # minval<- min(minval[minval!=0])
  # minval<- min(minval,1e-15)
  phi<- sapply(1:K,function(k)
    archmCopulaLink(copulafam=copulafam,param=alpha, p=S_k_D[[k]]))
  phiallK<- apply(phi,1,sum)
  
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
    
    
    Js<- sum(diff(c(SD.t2,0))*phi.t* phiprod.t)*(-1)^{dd+1+length(loc0)}
    
    
    
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
      # Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t)*(-1)^{dd+1}
      
      
      
      psi_d<- sapply(1:ny, function(ti)
        mean(exp(-xsim*phi1.t[ti])*exp(-xsim*phi00.tt[ti])*(-xsim)^dd,na.rm = TRUE))
      
      Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t)*(-1)^{dd+1}
      if(!is.na(Js)){
        if(Js>1){
          psi_d<- sapply(1:ny, function(ti)
            median(exp(-xsim*phi1.t[ti])*exp(-xsim*phi00.tt[ti])*(-xsim)^dd,na.rm = TRUE))
          
          Js<- sum(diff(c(SD.t2,0))*psi_d* phiprod.t)*(-1)^{dd+1}
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
      
      
      Js<- sum(diff(c(SD.t2,0))*psi_d)*(-1)
      
      
    
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
      
      
      Js<- sum(diff(c(SD.t2,0))*phi.t)*(-1)^{1+length(loc0)}
      
      
      
    }
  
  return(Js)
}


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

psialpha_dev<- function(d,copulafam,param){
  # psialpha_dev<- function(u,d,copulafam,param){
  # if(copulafam=="amh"){
  #   phi_inv<-  "(param-1)/(param - exp(y))"
  # }
  # if(copulafam=="clayton"){
  #   phi_inv<- "(param*y+1)^{-1/param}"
  # }
  # if(copulafam=="frank"){
  #   phi_inv<- "-1/param*log(1+(exp(-param)-1)*exp(y))"
  # }
  # if(copulafam=="gumbel"){
  #   phi_inv<- "exp(-y^(1/param))"
  # }
  # if(copulafam=="joe"){
  #   phi_inv<-"1- (1-exp(-y))^(1/param)"
  # }
  # calculus::derivative(f = phi_inv, var = c(y = u), 
  #                      order = d,params = list(param=param))
  
  
  if(copulafam=="amh"){
    phi_inv<- expression((param-1)/(param - exp(y)))
  }
  if(copulafam=="clayton"){
    phi_inv<- expression((param*y+1)^{-1/param})
  }
  if(copulafam=="frank"){
    phi_inv<- expression(-1/param*log(1+(exp(-param)-1)*exp(y)))
  }
  if(copulafam=="gumbel"){
    phi_inv<- expression(exp(-y^(1/param)))
  }
  if(copulafam=="joe"){
    phi_inv<-expression(1- (1-exp(-y))^(1/param))
  }
  
  df<-DD(phi_inv,"y",d)
  
  
  df
}
DD <- function(expr, name, order = 1,...) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name,...)
  else DD(D(expr, name,...), name, order - 1,...)
}

get_subsets<- function(vec){
  vec<- as.numeric(unique(vec))
  unlist(lapply(1:length(vec),   
                combn, 
                x = vec,
                simplify = FALSE), 
         recursive = FALSE)
  
}

Calcondsurv_D<- function(surv1,surv2,copulafam,copulaparam){
  if(length(surv1)>1&length(surv2)>1){
    out<- lapply(surv2, function(xx)
      dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =xx,mode = "2"))
  }else{
    out<- dev_Copula(copulafam = copulafam,param = copulaparam,p1 =surv1,p2 =surv2,mode = "2")
  }
  return(out)
  
} 


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
    
    sdev<- approxfun(x=x.new,y=sdev,method = "linear",
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
    
    
    sdev<- stepfun(x=x.new,y=c(0,diff(c(1,y.new))), right = FALSE,ties = mean )
  }
  
  
  
  return(sdev)
  
}