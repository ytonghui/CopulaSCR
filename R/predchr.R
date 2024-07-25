
predjsurv<- function(data,fit,t1.formula, t2.formula,t1, t2,method=c("np","sp")){
  Call <- match.call()
  method<- match.arg(method)
  
  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2'), 
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  # if(length(t1)!=length(t2)) stop("t1 and t2 shall have equal length")
  
  if(method=="np"){
    if (indx[1]==0&indx[2]==0) stop("a data argument is required")
    if (indx[3]==0&indx[2]==0) stop("a t1.formula argument is required")
    if (indx[4]==0&indx[2]==0) stop("a t2.formula argument is required")
    
    if(indx[2]!=0) {
      if(!class(fit)=="scrsurv"){stop("fit shall be 'scrsurv' class.")}
      t1.formula<- fit$t1.formula
      t2.formula<- fit$t2.formula
    }
    if(indx[1]!=0){
      t1data <- prepData(formula=t1.formula, data = data)
      t2data <- prepData(formula=t2.formula, data = data)
      
    }else{
      t1data<- fit$t1data
      t2data<- fit$t2data
      
    }
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
    return(list(js=js,t1=t1,t2=t2))
  }
    
  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if(class(fit)!="scrsurv"){stop("fit shall be 'scrsurv' class.")}
    
    surv1<- fit$t1.surv
    surv2<- fit$t2.surv
    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam
    
    s1.t1<- approx(x=c(0,surv1$time),y=c(1,surv1$surv),xout = t1,method = "constant",
                   yleft = 1,yright = min(surv1$surv),f = 0,na.rm = TRUE)$y
    
    s2.t2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = t2,method = "constant",
                   yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE)$y
    
    
    
    js<- lapply(1:length(t1),function(k1)
      Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.t2))
    names(js)<-  paste0("t1=",t1)
    return(list(js=js,t1=t1,t2=t2))
  }
}

predchr<- function(data,fit,t1.formula, t2.formula,tau, 
                   copulafam=c("clayton","frank","joe","gumbel","amh"),
                   t1, t2,method=c("np","sp")){
  
  Call <- match.call()
  method<- match.arg(method)
  copulafam<- match.arg(copulafam)
  
  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2', 'tau'), 
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  
  
  if(method=="np"){
    if (indx[1]==0&indx[2]==0) stop("a data argument is required")
    if (indx[3]==0&indx[2]==0) stop("a t1.formula argument is required")
    if (indx[4]==0&indx[2]==0) stop("a t2.formula argument is required")
    if (indx[7]==0&indx[2]==0) stop("a tau argument is required")
    
    if(indx[2]!=0) {
      if(!class(fit)=="scrsurv"){stop("fit shall be 'scrsurv' class.")}
      js<- predjsurv(fit = fit, t1=t1, t2=t2,method="np")
      copulafam<- fit$copulafam
      copulaparam<- fit$copulaparam
      
    }else{
      if(copulafam=="frank")copulaparam<- iTau(frankCopula(),tau)
      if(copulafam=="gumbel")copulaparam<- iTau(gumbelCopula(),tau)
      if(copulafam=="clayton")copulaparam<-iTau(claytonCopula(),tau)
      if(copulafam=="joe")copulaparam<-iTau(joeCopula(),tau)
      if(copulafam=="amh")copulaparam<-iTau(amhCopula(),tau)
      js<- predjsurv(data = data,t1.formula = t1.formula, t2.formula = t2.formula,
                     t1=t1, t2=t2,method="np")
    }
    
    
  }
  
  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if(class(fit)!="scrsurv"){stop("fit shall be 'scrsurv' class.")}
    
    
    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam
    js<- predjsurv(fit = fit, t1=t1, t2=t2,method="sp")
    
    
  }
  
  chr<- lapply(js$js,function(xx)
    npCHRfit(s=xx,param = copulaparam,copulafam = copulafam))
  
  return(list(chr=chr,t1=t1,t2=t2))
  
  
}



predcsurv<- function(data,fit,t1.formula, t2.formula,t1, t2, s2, 
                       t1equal =FALSE,method=c("np","sp")){
  Call <- match.call()
  method<- match.arg(method)
  
  indx <- match(c('data','fit','t1.formula', 't2.formula' ,'t1', 't2','s2'), 
                names(Call), nomatch=0)
  if (indx[5]==0) stop("a t1 argument is required")
  if (indx[6]==0) stop("a t2 argument is required")
  if (indx[7]==0) stop("a s2 argument is required")
  
  if(method=="np"){
    if (indx[1]==0&indx[2]==0) stop("a data argument is required")
    if (indx[3]==0&indx[2]==0) stop("a t1.formula argument is required")
    if (indx[4]==0&indx[2]==0) stop("a t2.formula argument is required")
    
    if(indx[2]!=0) {
      if(!class(fit)=="scrsurv"){stop("fit shall be 'scrsurv' class.")}
      t1.formula<- fit$t1.formula
      t2.formula<- fit$t2.formula
      copulaparam<- fit$copulaparam
    }
    if(indx[1]!=0){
      t1data <- prepData(formula=t1.formula, data = data)
      t2data <- prepData(formula=t2.formula, data = data)
      
    }else{
      t1data<- fit$t1data
      t2data<- fit$t2data
      
    }
    
    
    
    if(isTRUE(t1equal))stop("t1equal shall be FALSE for method 'np'.")
    
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
    return(list(condsurv = condsurv,t1=t1,t2=t2,s2=s2))
  }
  if(method=="sp"){
    if (indx[2]==0) stop("a fit argument is required")
    if(class(fit)!="scrsurv"){stop("fit shall be 'scrsurv' class.")}
    
    surv1<- fit$t1.surv
    surv2<- fit$t2.surv
    copulafam<- fit$copulafam
    copulaparam<- fit$copulaparam
    
    s1.t1<- approx(x=c(0,surv1$time),y=c(1,surv1$surv),xout = t1,method = "constant",
                   yleft = 1,yright = min(surv1$surv),f = 0,na.rm = TRUE,ties = mean)$y
    
    s2.t2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = t2,method = "constant",
                   yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE,ties = mean)$y
    
    s2.s2<- approx(x=c(0,surv2$time),y=c(1,surv2$surv),xout = s2,method = "constant",
                   yleft = 1,yright = min(surv2$surv),f = 0,na.rm = TRUE,ties = mean)$y
    if(!isTRUE(t1equal)){
      js.t<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.t2))
      
      
      js.s<- lapply(1:length(t1),function(k1)
        Copulafn(copulafam = copulafam,param = copulaparam,p1 = s1.t1[k1],p2 = s2.s2))
    }
    if(isTRUE(t1equal)){
      # js.t<-  lapply(1:length(t1),function(k)
      #   archmCopulaLink_dev(copulafam=copulafam,param=copulaparam, js=js.t[[k]]))
      js.t<-  lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param = copulaparam,p1 =s1.t1[k1],p2 =s2.t2,mode = "1"))
        
      # js.s<-  lapply(1:length(t1),function(k)
      #   archmCopulaLink_dev(copulafam=copulafam,param=copulaparam, js=js.s[[k]]))
      
      js.s<- lapply(1:length(t1),function(k1)
        dev_Copula(copulafam = copulafam,param = copulaparam,p1 =s1.t1[k1],p2 =s2.s2,mode = "1"))
      
    }
    names(js.t)<- names(js.s)<- paste0("t1=",t1)
    
    condsurv<- lapply(1:length(t1), function(k)tcrossprod(js.t[[k]],1/js.s[[k]]))
    condsurv<- lapply(condsurv,t2s2matNA,t2=t2,s2=s2)
    
    condsurv<- lapply(1:length(t1), function(k)t1s2matNA(x=condsurv[[k]],tt1=t1[k],s2))
    names(condsurv)<- paste0("t1=",t1)
    return(list(condsurv = condsurv,t1=t1,t2=t2,s2=s2))
  }
  
}


t2s2matNA<- function(x,t2,s2){
  
  colnames(x)<- paste0("s2=",s2)
  rownames(x)<- paste0("t2=",t2)
  
  
  tsind<- sapply(s2,function(xx)xx>t2)
  x[tsind]<- NA
  return(x)
}


t1s2matNA<- function(x,tt1,s2){
  x[,tt1>s2]<- NA
  return(x)
}


