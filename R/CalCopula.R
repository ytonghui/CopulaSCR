#' @title Dependence Measures for Bivariate Copulas
#' @description  Compute Kendall's Tau of an Archimedean copula and
#' copula parameter given the value of Kendall's Tau.
#' @aliases Caltau
#' @aliases Calitau
#'
#'
#' @usage Caltau(copulafam, copulaparam)
#' @usage Calitau(copulafam, tau)
#' @param copulafam a character string specifying the family of an Archimedean copula.
#' Currently supported families are "frank", "clayton",  "amh", "gumbel", and "joe".
#' @param copulaparam number (numeric) specifying the copula parameter.
#' @param tau Kendall's tau of an Archimedean copula.
#'
#'
#' @examples Caltau(copulafam = "frank", copulaparam = 2) # output: tau
#' Calitau(copulafam = "frank", tau = 0.5) # output: copulaparam
#'
#' @importFrom copula tau iTau
## #' @import acopula
#' @seealso \code{\link[copula:tau]{tau}}
#' @export Caltau
#' @export Calitau

Caltau<- function(copulafam,copulaparam){
  if(copulafam=="frank")tau<- copula::tau(copula::frankCopula(copulaparam))
  if(copulafam=="gumbel")tau<- copula::tau(copula::gumbelCopula(copulaparam))
  if(copulafam=="clayton")tau<-copula::tau(copula::claytonCopula(copulaparam))
  if(copulafam=="joe")tau<- copula::tau(copula::joeCopula(copulaparam))
  if(copulafam=="amh")tau<- copula::tau(copula::amhCopula(copulaparam))

  return(tau)
}

Calitau<- function(copulafam,tau){
  if(copulafam=="frank")alpha<- copula::iTau(copula::frankCopula(),tau)
  if(copulafam=="gumbel")alpha<- copula::iTau(copula::gumbelCopula(),tau)
  if(copulafam=="clayton")alpha<-copula::iTau(copula::claytonCopula(),tau)
  if(copulafam=="joe")alpha<- copula::iTau(copula::joeCopula(),tau)
  if(copulafam=="amh")alpha<- copula::iTau(copula::amhCopula(),tau)
  return(alpha)
}

#' @noRd
archmCopulaLink<- function(copulafam,param, p){
  # archmCopula includes  clayton, frank, amh, gumbel, joe
  if(copulafam=="amh"){
    phi<- log(param +(1-param)/p)
  }
  if(copulafam=="clayton"){
    phi<- p^{-param}-1
  }
  if(copulafam=="frank"){
    phi<- -log((exp(-param*p)-1)/(exp(-param)-1))
  }
  if(copulafam=="gumbel"){
    phi<- (-log(p))^param
  }
  if(copulafam=="joe"){
    phi<-  -log( 1- (1-p)^param)
  }
  return(phi)
}

#' @noRd
archmCopulaLink_inv<- function(copulafam,param, y){
  # archmCopula includes  clayton, frank, amh, gumbel, joe
  if(copulafam=="amh"){
    phi_inv<-  (param-1)/(param - exp(y))
  }
  if(copulafam=="clayton"){
    phi_inv<- (y+1)^{-1/param}
  }
  if(copulafam=="frank"){
    phi_inv<- -1/param*log(1+(exp(-param)-1)*exp(-y))
  }
  if(copulafam=="gumbel"){
    phi_inv<- exp(-y^(1/param))
  }
  if(copulafam=="joe"){
    phi_inv<-1- (1-exp(-y))^(1/param)
  }
  return(phi_inv)

}

#' @noRd
archmCopulaLink_dev<- function(copulafam,param, js){
  # archmCopula includes  clayton, frank, amh, gumbel, joe
  if(copulafam=="amh"){
    phi_dev<-  (param-1)/(param*js^2+(1-param)*js)
  }
  if(copulafam=="clayton"){
    phi_dev<- - js^{-param-1}*param
  }
  if(copulafam=="frank"){
    phi_dev<- param* exp(-param*js)/(exp(-param*js)-1)
  }
  if(copulafam=="gumbel"){
    phi_dev<- -param*((-log(js))^{param-1})/js
  }
  if(copulafam=="joe"){
    phi_dev<- -param*((1-js)^{param-1})/(1-(1-js)^param)
  }
  return(phi_dev)

}

#' @noRd
Copulafn<- function(copulafam,param, p1,p2){
  # archmCopula includes  clayton, frank, amh, gumbel, joe

  a<- switch(copulafam,
             amh = {
               if(param>=-1 &param<=1) p1*p2/(1-param*(1-p1)*(1-p2))
               else 0

               # pCopula(c(p1, p2), amhCopula(param = param,dim=2,use.indepC="FALSE"))
             },
             clayton = {
               if(param>= -1)(p1^{-param}+p2^{-param}-1)^(-1/param)
               else 0
               # pCopula(c(p1, p2), claytonCopula(param = param,dim=2,use.indepC="FALSE"))
             },
             frank = {
               if(param!=0) -log(1+(exp(-param*p1)-1)*(exp(-param*p2)-1)/(exp(-param)-1))/param
               else 0
               # pCopula(c(p1, p2), frankCopula(param = param,dim=2,use.indepC="FALSE"))
             },
             gumbel = {
               if(param>=1) exp(-((-log(p1))^param+(-log(p2))^param)^{1/param})
               else 0

               # pCopula(c(p1, p2), gumbelCopula(param = param,dim=2,use.indepC="FALSE"))
             },
             joe = {
               if(param>=1) 1-((1-p1)^param+(1-p2)^param-((1-p1)^param)*((1-p2)^param))^(1/param)
               else 0
               # pCopula(c(p1, p2), joeCopula(param = param,dim=2,use.indepC="FALSE"))
             },
  )
  a[is.na(a)]<- 0
  return(a)
}

#' @noRd
dev_Copula<- function(copulafam, param, p1,p2, mode = "1"){

  if(mode=="1"&copulafam=="amh"){
    if(param>=-1 &param<=1) {
      a<- (1-param*(1-p1)*(1-p2))
      a<- (p2*a- p1*p2*param*(1-p2))/a^2
    }
    else{
      a<- NA
    }
    return(a)
  }else  if(mode=="2"&copulafam=="amh"){
    if(param>=-1 &param<=1) {
      a<- (1-param*(1-p1)*(1-p2))
      a<- (p1*a- p1*p2*param*(1-p1))/a^2
    }
    else{
      a<- NA
    }
    return(a)
  }else  if(mode=="12"&copulafam=="amh"){
    if(param>=-1 &param<=1) {
      a<- (1-param*(1-p1)*(1-p2))
      a<- 1-param+2*param*p2/a^2 - 2*param*p2*(1-param+param*p2)*(1-p1)/a^3
    }else{
      a<- NA
    }
    return(a)
  }else  if(mode=="1"&copulafam=="clayton"){

    if(param>= -1) {
      out<- (p1^{-param}+p2^{-param}-1)^(-1/param-1)
      out<- out*p1^{-param-1}
      out[is.na(out)]<- 0
    }  else {
      out<- NA
    }
    return(out)


  }else  if(mode=="2"&copulafam=="clayton"){

    if(param>= -1) {
      out<- (p1^{-param}+p2^{-param}-1)^(-1/param-1)
      out<- out*p2^{-param-1}
      out[is.na(out)]<- 0
    }  else {
      out<- NA
    }
    return(out)

  }else  if(mode=="12"&copulafam=="clayton"){
    if(param>= -1) {
      out<- (p1^{-param}+p2^{-param}-1)^(-1/param-2)
      out<- (1+param)*out*p1^{-param-1}*p2^{-param-1}
    }else{
      out<- NA
    }
    return(out)
  }else  if(mode=="1"&copulafam=="frank"){

    if(param!=0) {
      out1<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out<- (exp(-param*p2)-1)*exp(-param*p1)/((exp(-param)-1)+out1)

    } else {
      out<- NA
    }
    return(out)
  }else  if(mode=="2"&copulafam=="frank"){

    if(param!=0) {
      out1<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out<- (exp(-param*p1)-1)*exp(-param*p2)/((exp(-param)-1)+out1)

    }  else {out<- NA}
    return(out)

  }else  if(mode=="11"&copulafam=="frank"){

    if(param!=0) {
      out0<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out1<- out0+(exp(-param)-1)
      out<- param*(exp(-param*p2)-1)*exp(-param*p1)*(exp(-param*p2)-exp(-param))/(out1^2)

    }    else {out<-NA}

    return(out)

  }else  if(mode=="12"&copulafam=="frank"){

    if(param!=0) {
      out0<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out1<- out0+(exp(-param)-1)
      out<- (1-exp(-param))*param*exp(-param*p1)*exp(-param*p2)/(out1^2)

    } else{
      out<- NA
    }
    return(out)
  }else  if(mode=="22"&copulafam=="frank"){

    if(param!=0) {
      out0<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out1<- out0+(exp(-param)-1)
      out<- param*(exp(-param*p1)-1)*exp(-param*p2)*(exp(-param*p1)-exp(-param))/(out1^2)

    } else {out<- NA}

    return(out)


  }else  if(mode=="121"&copulafam=="frank"){

    if(param!=0) {
      out0<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out1<- out0+(exp(-param)-1)
      out2<- (1-exp(-param))*param*exp(-param*p1)*exp(-param*p2)
      out<- out2*(-param)*(out1^2)-2*out2*out1*(exp(-param*p2)-1)*exp(-param*p1)*(-param)
      out<- out/(out1^4)

    }else {
      out<- NA
    }
    return(out)

  }else  if(mode=="122"&copulafam=="frank"){


    if(param!=0) {
      out0<- (exp(-param*p1)-1)*(exp(-param*p2)-1)
      out1<- out0+(exp(-param)-1)
      out2<- (1-exp(-param))*param*exp(-param*p1)*exp(-param*p2)
      out<- out2*(-param)*(out1^2)-2*out2*out1*(exp(-param*p1)-1)*exp(-param*p2)*(-param)
      out<- out/(out1^4)

    } else{out<-  NA
    }

    return(out)
  }else if(mode=="1"&copulafam=="gumbel"){

    if(param>=1) {
      out<- exp(-((-log(p1))^param+(-log(p2))^param)^{1/param})*
        ((-log(p1))^param+(-log(p2))^param)^{1/param-1}*(-log(p1))^{param-1}/p1
      # if(is.na(out)&p1==1&p2==1)out<- 0
    }  else {
      out<- NA
    }
    return(out)

  }else  if(mode=="2"&copulafam=="gumbel"){

    if(param>=1) {
      out<- exp(-((-log(p1))^param+(-log(p2))^param)^{1/param})*
        ((-log(p1))^param+(-log(p2))^param)^{1/param-1}*(-log(p2))^{param-1}/p2
      # if(is.na(out)&p1==1&p2==1)out<- 0
    }  else {
      out<- NA
    }
    return(out)


  }else  if(mode=="12"&copulafam=="gumbel"){

    if(param>=1) {
      out<- exp(-((-log(p1))^param+(-log(p2))^param)^{1/param})*
        { ((-log(p1))^param+(-log(p2))^param)^{2/param-2}*(-log(p1))^{param-1}/p1*
            (-log(p2))^{param-1}/p2 +
            ((-log(p1))^param+(-log(p2))^param)^{1/param-2}*(-log(p1))^{param-1}/p1*
            (-log(p2))^{param-1}/p2*(param-1)}
      # if(is.na(out)&p1==1&p2==1)out<- 0
    }  else {
      out<- NA
    }
    return(out)

  }else if(mode=="1"&copulafam=="joe"){

    if(param>=1) {
      a<- ((1-p1)^param+(1-p2)^param-((1-p1)^param)*((1-p2)^param))^(1/param-1)
      a<- a*(((1-p1)^(param-1))-((1-p2)^param)*((1-p1)^(param-1)))
    }    else {a<-NA}

    return(a)

  }else  if(mode=="2"&copulafam=="joe"){

    if(param>=1) {
      a<- ((1-p1)^param+(1-p2)^param-((1-p1)^param)*((1-p2)^param))^(1/param-1)
      a<- a*(((1-p2)^(param-1))-((1-p1)^param)*((1-p2)^(param-1)))
    } else {
      a<- NA
    }
    return(a)

  }else  if(mode=="12"&copulafam=="joe"){

    if(param>=1) {
      a<- (1-p1)^param+(1-p2)^param-((1-p1)^param)*((1-p2)^param)
      a<- (param-1)*a^(1/param-2)*((1-p1)^(param-1))*((1-p2)^(param-1))*
        (1-(1-p1)^param)*(1-(1-p2)^param)+
        a^(1/param-1)*param*((1-p1)^(param-1))*((1-p2)^(param-1))

    }    else {a<-NA}
    return(a)
  }else{
    if(mode=="1") devmode <- c(1,0)
    if(mode=="2") devmode <- c(0,1)
    if(mode=="11") devmode <- c(2,0)
    if(mode=="12") devmode <- c(1,1)
    if(mode=="22") devmode <- c(0,2)
    if(mode=="121") devmode <- c(2,1)
    if(mode=="122") devmode <- c(1,2)
    out<- acopula::nderive(fun = function(x) {
      u<- x[1]
      v<- x[2]
      Copulafn(copulafam =copulafam ,param = param, p1=u,p2=v)
    }, point = c(p1,p2),order = devmode, difference = 1e-06)
    return(out)
  }


}






copulabd<- function(copulafam){
  lb<- switch(copulafam,
              amh = -1,
              clayton = -1,
              frank = -1000,
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
