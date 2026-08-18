#' Dependence Measures for Bivariate Copulas
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
#' @return A numeric value or numeric vector representing Kendall's tau for the
#'   specified copula family and association parameter.
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
archmCopulaLink<- function(copulafam, param, p) {
  archmCopulaLink_cpp(copulafam, param, p)
}

#' @noRd
archmCopulaLink_dev<- function(copulafam, param, js) {
  archmCopulaLink_dev_cpp(copulafam, param, js)
}

#' @noRd
archmCopulaLink_inv<- function(copulafam, param, y) {
  archmCopulaLink_inv_cpp(copulafam, param, y)
}

#' @noRd
dev_Copula<- function(copulafam, param, p1, p2, mode = "1") {
  dev_Copula_cpp(copulafam, param, p1, p2, mode)
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
