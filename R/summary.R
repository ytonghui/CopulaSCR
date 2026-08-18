#' Summarize an 'scrassonp' Object
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
    stratas <- object$stratas
    strata.names <- as.character(stratas[[1]])
  }else{
    stratas <- NULL
    strata.names <- NULL
  }

  if(isTRUE(object$se)){
    paramses<- object$param.se
    tauses<- object$tau.se
    if(length(taus)>1) {
      params.lower<- apply(object$param.boot,2,quantile,probs = (1-conf.int)/2)
      params.upper<- apply(object$param.boot,2,quantile,probs = 1- (1-conf.int)/2)

      taus.lower<- apply(object$tau.boot,2,quantile,probs = (1-conf.int)/2)
      taus.upper<- apply(object$tau.boot,2,quantile,probs = 1- (1-conf.int)/2)
      taus.upper <- pmin(taus.upper, 1)

      # paramout<- as.matrix(data.frame(param = params,param.StdErr=paramses,
      #                                 param.lower =params.lower, param.upper =params.upper))
      #
      # tauout<- as.matrix(data.frame(tau = taus,tau.StdErr = tauses,
      #                               tau.lower =taus.lower, tau.upper =taus.upper))

      paramout<- cbind(params,paramses,params.lower, params.upper)

      tauout<- cbind(taus,tauses,taus.lower,taus.upper)


      rownames(paramout) <- rownames(tauout) <- strata.names
      colnames(paramout)<- colnames(tauout)<-
        c("Est","StdErr",paste0(conf.int*100,"%lower"), paste0(conf.int*100,"%upper"))

    }else{
      params.lower<- quantile(object$param.boot,probs = (1-conf.int)/2)
      params.upper<- quantile(object$param.boot,probs = 1- (1-conf.int)/2)
      taus.lower<- quantile(object$tau.boot,probs = (1-conf.int)/2)
      taus.upper<- quantile(object$tau.boot,probs = 1- (1-conf.int)/2)
      taus.upper<- pmin(taus.upper, 1)
      


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
      names(params) <- names(taus) <- strata.names
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


#' Summarize an 'mscr' Object
#'
#' Summarizing the association analysis for semi-competing risks data
#' with multiple intermediate event times.
#'
#' @param object An object of class \code{mscr} derived with
#' \code{\link{mscr}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{summary.mscr} containing the
#' estimated copula parameters and Kendall's tau values.
#'
#' @seealso \code{\link{mscr}}, \code{\link{print.mscr}}
#'
#' @keywords survival
#' @export
summary.mscr <- function(object, ...) {
  
  theta <- object$theta
  tau.theta <- object$tau.theta
  
  if (is.null(names(theta))) {
    names(theta) <- paste0("T1.", seq_along(theta))
  }
  
  if (is.null(names(tau.theta))) {
    names(tau.theta) <- names(theta)
  }
  
  sumout <- list(
    Call = object$call,
    copulafam = object$copulafam,
    theta = theta,
    alpha = object$alpha,
    tau.theta = tau.theta,
    tau.alpha = object$tau.alpha
  )
  
  class(sumout) <- "summary.mscr"
  sumout
}


#' @export
print.summary.mscr <- function(x, ...) {
  
  cat("\n")
  cat("Call: ")
  print(x$Call)
  cat("\n")
  
  cat("Association analysis of semi-competing risks data with multiple intermediate events.\n")
  
  cat("Estimator of copula parameters in the", toupper(x$copulafam),
      "copula structure\n")
  
  cat("theta:\n")
  print(round(x$theta, 3))
  
  cat("alpha:\n")
  print(round(x$alpha, 3))
  
  cat("----\n")
  
  cat("Estimator of Kendall's tau in the", toupper(x$copulafam),
      "copula structure\n")
  
  cat("tau.theta:\n")
  print(round(x$tau.theta, 3))
  
  cat("tau.alpha:\n")
  print(round(x$tau.alpha, 3))
  
  invisible(x)
}


#' Summarize an 'scrsurv' Object
#'
#' Summarizing the estimated marginal survival functions for the
#' nonterminal and terminal event times.
#'
#' @param object An object of class \code{scrsurv} derived with
#' \code{\link{scrsurv}}.
#' @param censored Logical. If \code{TRUE}, censoring times are included
#' in the output. Otherwise, only event times are included.
#' @param scale Numeric value used to rescale the survival times.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{summary.scrsurv} containing summaries
#' of the marginal survival functions for the nonterminal and terminal
#' event times.
#'
#' @seealso \code{\link{scrsurv}}, \code{\link{plot.scrsurv}}
#'
#' @keywords survival
#' @export
summary.scrsurv <- function(object, censored = FALSE, scale = 1, ...) {
  
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("'scale' must be a positive numeric value.")
  }
  
  make_summary <- function(x) {
    
    keep <- if (isTRUE(censored)) {
      rep(TRUE, length(x$time))
    } else {
      x$n.event > 0
    }
    
    out <- data.frame(
      time = x$time[keep] / scale,
      n.risk = x$n.risk[keep],
      n.event = x$n.event[keep],
      n.censor = x$n.censor[keep],
      surv = x$surv[keep]
    )
    
    if (!is.null(x$std.err)) {
      out$std.err <- x$std.err[keep]
    }
    
    if (!is.null(x$lower)) {
      out$lower <- x$lower[keep]
    }
    
    if (!is.null(x$upper)) {
      out$upper <- x$upper[keep]
    }
    
    out
  }
  
  sumout <- list(
    Call = object$Call,
    t1.surv = make_summary(object$t1.surv),
    t2.surv = make_summary(object$t2.surv),
    na.action = object$na.action
  )
  
  class(sumout) <- "summary.scrsurv"
  sumout
}


#' @export
print.summary.scrsurv <- function(x, ...) {
  
  cat("\n")
  cat("Call: ")
  print(x$Call)
  cat("\n")
  
  cat("Nonterminal event survival curve:\n")
  print(x$t1.surv, row.names = FALSE)
  
  cat("\n")
  cat("Terminal event survival curve:\n")
  print(x$t2.surv, row.names = FALSE)
  
  if (!is.null(x$na.action)) {
    cat("\n",
        length(x$na.action),
        ifelse(length(x$na.action) == 1,
               " observation",
               " observations"),
        " deleted due to missing values.\n",
        sep = "")
  }
  
  invisible(x)
}