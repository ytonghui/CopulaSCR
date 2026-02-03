#' @noRd
predictMSCR<- function(fit, t1obs,times,type=c("surv", "rrms", "rmst", "qrl","qst"),
                       tu,tau=0.5,nsim=1000,smooth= FALSE,jumpsize=FALSE,exact=TRUE){
  #times	：Vector of times at which to return the estimated survival rates.
  Call <- match.call()
  indx <- match(c('fit', 't1obs','times','tu'), names(Call), nomatch=0)
  type<- match.arg(type)

  if (indx[1]==0) stop("a fit argument of class 'mscr' is required")
  if (indx[2]==0) stop("a t1obs argument is required")


  if(!inherits(fit,"mscr")){
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


#' @noRd
predbsMSCR<- function(fit, t1obs,times,type=c("surv","rrms", "rmst", "qrl","qst"),
                      tu,tau=0.5,cause=0,maxobs=TRUE){
  #times	：Vector of times at which to return the estimated survival rates.
  Call <- match.call()
  indx <- match(c('fit', 't1obs','times','tu'), names(Call), nomatch=0)
  type<- match.arg(type)

  if (indx[1]==0) stop("a fit argument of class 'mscr' is required")
  if (indx[2]==0) stop("a t1obs argument is required")


  if(!inherits(fit,"mscr")){
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



###### auxillary functions ###########

#' @noRd
dypredfit<- function(fit,t1obs,type,times,nsim,tu,tau,smooth,jumpsize,exact=FALSE){
  ## type = "surv": survival function
  ## type = "rrms": restricted residual mean survival
  ## type = "rmst": restricted residual mean survival+ landmark time rrms+max(t1obs)
  ## type = "qrl": quantile residual lifetime
  ## type = "qst": quantile survival time (qrl+ landmark time)

  if(all(is.na(t1obs))){
    return(bspredfit(fit = fit,t1obs = t1obs,type = type,times = times,
                     tu = tu,tau = tau,cause = 0))
  }
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
    cs<- predcsurv(fit=fit$mar.fits[[t1pos]],t1=t1obs, t2=times, s2=t1obs,
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

#' @noRd
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
      cs<- predcsurv(fit=fit$mar.fits[[t1pos]],t1=t1obs2, t2=times, s2=tmax,
                     t1equal =TRUE,method="sp")$condsurv[[1]][,1]
    }else{
      cs<- predcsurv(fit=fit$mar.fits[[t1pos]],t1=t1obs2, t2=times, s2=t1obs2,
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



#' @noRd
dypred.Q<- function(fit,t1obs,t1pos,nsim=500,smooth,jumpsize,exact=FALSE){
  jps<- knots(fit$S_D)
  SD.jp<- sapply(jps, fit$S_D)
  lpos<- length(t1obs)
  njp<- length(jps)

  Sk.T<- sapply(1:lpos, function(k)fit$S_k[[t1pos[k]]](t1obs[k]))

  # if(!isTRUE(smooth)){
  # Sdev_k<- lapply(1:lpos,function(k)
  #   stepfun(x=fit$mar.fits[[t1pos[k]]]$t1.surv$time,
  #           y=c(0,diff(c(1,fit$mar.fits[[t1pos[k]]]$t1.surv$surv))),
  #           right = FALSE))

  # }else{
  # Sdev_k<- lapply(1:lpos,function(k)
  #   data.frame(x=fit$mar.fits[[t1pos[k]]]$t1.surv$time,
  #              y=diff(c(1,fit$mar.fits[[t1pos[k]]]$t1.surv$surv))))
  # Sdev_k<- lapply(1:lpos,function(k)Sdev_k[[k]][Sdev_k[[k]]$y!=0,])
  # Sdev_k<- lapply(1:lpos,function(k)
  #   stepfun(x=Sdev_k[[k]]$x,y=c(0,Sdev_k[[k]]$y), right = FALSE))
  # }

  Sdev_k<- lapply(1:lpos,function(k)
    survnumdiff(fit$mar.fits[[t1pos[k]]]$t1.surv$time,
                fit$mar.fits[[t1pos[k]]]$t1.surv$surv,
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







