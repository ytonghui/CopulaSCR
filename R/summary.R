#' Summary method for scrassonp objects.
#'
#' Summarizing the result of the association analysis in the semi-competing risks data.
#' Confidence intervals are displayed when they are part of the fitted
#' object.
#'
#' @param object An object with class `scrassonp' derived with
#' \code{\link{scrassonp}}
#' @param conf.int confidence level
#' @param ... Further arguments that are passed to the print
#' function.
#' @return A data.frame with the relevant information.
#' @seealso \code{\link{scrassonp}}
#'
#' @keywords survival
#' @export
summary.scrassonp <- function(object,conf.int=TRUE,...) {

  if (isTRUE(conf.int)|isTRUE(object$se)) conf.int <- 0.95

  if (is.numeric(conf.int)) {
    if (!(0 < conf.int && conf.int < 1))
      conf.int <- 0.95
  }



  params<- object$copulaparam
  taus<- object$tau

  if(length(taus)>1) {
    stratas<- object$t1data$stratas
  }else{
    stratas<- NULL
  }

  if(isTRUE(object$se)){
    paramses<- object$param.se
    tauses<- object$tau.se
    if(length(taus)>1) {
      params.lower<- apply(object$param.boot,2,quantile,probs = (1-conf.int)/2)
      params.upper<- apply(object$param.boot,2,quantile,probs = 1- (1-conf.int)/2)

      taus.lower<- apply(object$tau.boot,2,quantile,probs = (1-conf.int)/2)
      taus.upper<- apply(object$tau.boot,2,quantile,probs = 1- (1-conf.int)/2)

      # paramout<- as.matrix(data.frame(param = params,param.StdErr=paramses,
      #                                 param.lower =params.lower, param.upper =params.upper))
      #
      # tauout<- as.matrix(data.frame(tau = taus,tau.StdErr = tauses,
      #                               tau.lower =taus.lower, tau.upper =taus.upper))

      paramout<- cbind(params,paramses,params.lower, params.upper)

      tauout<- cbind(taus,tauses,taus.lower,taus.upper)


      rownames(paramout)<- rownames(tauout)<- paste0("strata", 1:length(taus))
      colnames(paramout)<- colnames(tauout)<-
        c("Est","StdErr",paste0(conf.int*100,"%lower"), paste0(conf.int*100,"%upper"))

    }else{
      params.lower<- quantile(object$param.boot,probs = (1-conf.int)/2)
      params.upper<- quantile(object$param.boot,probs = 1- (1-conf.int)/2)
      taus.lower<- quantile(object$tau.boot,probs = (1-conf.int)/2)
      taus.upper<- quantile(object$tau.boot,probs = 1- (1-conf.int)/2)


      paramout<- c("param" = params,"param.StdErr"=paramses,
                   "param.lower" =params.lower, "param.upper" =params.upper)


      tauout<- c("tau" = taus,"tau.StdErr" = tauses,
                 "tau.lower" =taus.lower, "tau.upper" =taus.upper)
      names(paramout)<-names(tauout)<-
        c("Est","StdErr",paste0(conf.int*100,"%lower"), paste0(conf.int*100,"%upper"))
    }




    sumout<- list(Call =object$Call,copulafam = object$copulafam,
                  nstrata =length(taus),stratas=stratas,
                  param = paramout,tau = tauout)
  }else{
    if(length(taus)>1) {
      names(params)<- names(taus)<- paste0("strata", 1:length(taus))
    }
    sumout<- list(Call =object$Call,copulafam = object$copulafam,
                  nstrata =length(taus),stratas=stratas,
                  param = params,tau = taus)
  }

  class(sumout) <- "summary.scrassonp"
  sumout
}

#' @export
print.summary.scrassonp <- function(x,
                                    ...) {
  x<- x
  cat("\n")
  cat("Call: ")
  print(x$Call)
  cat("\n")


  if(x$nstrata==1){
    cat("Association analysis of semi-competing risks data.\n")
  }
  if(x$nstrata>1){
    cat("Stratified association analysis of semi-competing risks data.\n")
  }

  paramout<- x$param
  tauout<- x$tau
  cat("Estimator of copula parameter in the", toupper(x$copulafam),
      "copula structure \n")


  print(round(paramout,3))
  cat("----\n")
  cat("Estimator of Kendall's tau in the", toupper(x$copulafam),
      "copula structure \n")

  print(round(tauout,3))
  invisible(x)

}
