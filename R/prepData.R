#' @aliases prepData0
#' @noRd
#'
prepData0<- function (t1.formula,t2.formula, data){
  varnames<- union(all.vars(t1.formula),all.vars(t2.formula))

  data2<- na.omit(data[,varnames],varnames)
  na.action<- attr(data2, "na.action")
  return(list(data=data2,na.action=na.action))
}

#' @aliases prepData
#' @noRd
#'
prepData<- function (formula, data, na.action=NULL,
                     bandwidth = NULL, discrete.level = 10,xlimits=10)
{
  #discrete.level	:  Numeric covariates are treated as factors when their number
  #of unique values exceeds not discrete.level. Otherwise the product limit
  #method is applied, in overlapping neighborhoods according to the bandwidth.
  #x: 	  logical value: if TRUE, the full covariate matrix with is returned in
  #component model.matrix. The reduced matrix contains unique rows of the full
  #covariate matrix and is always returned in component X.

  # call <- match.call()

  if (!(is.matrix(data) | is.data.frame(data))) {
    stop("data object must be a matrix or a data.frame")
  }

  EHF <- prodlim::EventHistory.frame(formula = formula, data = data,
                                     unspecialsDesign = FALSE,
                                     specials = c("Strata", "strata", "factor", "NN", "cluster"),
                                     stripSpecials = c("strata", "cluster", "NN"),
                                     stripAlias = list(strata = c("Strata", "factor")),
                                     stripArguments = list(strata = NULL, NN = NULL, cluster = NULL),
                                     specialsDesign = FALSE, check.formula = TRUE)


  response <- EHF$event.history
  cens.type <- attr(response, "cens.type")
  if(cens.type!= "rightCensored"){
    stop("This algorithm can only accommodate right censored data.")
  }

  covariates <- EHF[-1]
  if (length(covariates$cluster) > 0) {
    stop("Cluster-correlated observations not yet handled.")
  }


  event.history <- EHF$event.history
  event.time.order <- 1:nrow(event.history)#order(event.history[, "time"], -event.history[, "status"])
  strata.pos <- match(c("strata", "factor", "Strata"), names(covariates),
                      nomatch = 0)
  if (sum(strata.pos) > 0) {
    strata <- do.call("cbind", covariates[strata.pos])
  }  else strata <- NULL
  NN <- covariates[["NN"]]
  xlevels <- attr(strata, "levels")
  rest <- covariates$design
  xlevels <- c(attr(strata, "levels"), attr(rest, "levels"))


  # cotype = 1: Surv(time,status)~ 1
  # cotype = 2: Surv(time,status)~ strata(sex)
  # cotype = 3: Surv(time,status)~ age
  # cotype = 4: Surv(time,status)~ strata(sex)+ age
  if ((is.null(NN) + is.null(strata) + is.null(rest)) == 3) {
    cotype <- 1
  }  else {
    unspecified <- NULL
    if (!is.null(rest)) {
      discrete.p <- sapply(colnames(rest), function(u) {
        x <- rest[, u, drop = TRUE]
        !is.numeric(x) || !length(unique(x)) > discrete.level
      })
      if (any(!discrete.p)) {
        NN <- if (is.null(NN))
          rest[, !discrete.p, drop = FALSE]
        else cbind(NN, rest[, !discrete.p, drop = FALSE])
      }
      if (any(discrete.p)) {
        strata <- if (is.null(strata)) {
          rest[, discrete.p, drop = FALSE]
        } else {
          cbind(strata, rest[, discrete.p, drop = FALSE])
        }
      }
    }
    if (NCOL(NN) > 1)
      stop(paste("Currently we can not compute neighborhoods in",
                 length(colnames(NN)), "continuous dimensions."))

    cotype <- 1 + (!is.null(strata)) * 1 + (!is.null(NN)) *
      2
  }

  if(cotype%in%c(3,4)){stop("Only one discrete variable is allowed.")}
  if(length(colnames(strata))>1){stop("Only one discrete variable is allowed.")}

  if (any(found <- (match(colnames(strata), names(xlevels),
                          nomatch = 0) == 0))) {
    uniquelevels <- lapply(colnames(strata)[found], function(x) {
      unique(strata[, x])
    })
    names(uniquelevels) <- colnames(strata)[found]
    xlevels <- c(xlevels, uniquelevels)
  }

  if (cotype %in% c(2, 4)) {
    S <- interaction(data.frame(strata), sep = ":", drop = TRUE)
    NS <- length(unique(S))
    Sfactor <- factor(S, levels = levels(S), labels = 1:NS)
    # sorted <- order(Sfactor, response[, "time"], -response[,  "status"])
    sorted <- order(Sfactor)
    Sfactor <- Sfactor[sorted]
  }  else {
    sorted <- event.time.order
  }
  response <- response[sorted, ]  ## sorted by treatment
  # if (is.null(weights)) {
  #   weighted <- 0
  #   weights <- NULL
  # }  else {
  #   weighted <- 1
  #   if (length(weights) != NROW(response))
  #     stop(paste("The length of caseweights is: ", length(weights),
  #                "\nthis is not the same as the number of subjects\n
  #                with no missing values, which is ",
  #                NROW(response), sep = ""))
  #   weights <- weights[sorted]
  # }
  if (cotype %in% c(3, 4)) {
    Z <- NN[sorted, , drop = TRUE]
    if (cotype == 3) {
      nbh <- prodlim::neighborhood(Z, bandwidth = bandwidth)
      nbh.list <- list(nbh)
      bandwidth <- nbh$bandwidth
      neighbors <- nbh$neighbors
    }    else {
      nbh.list <- lapply(split(Z, Sfactor), prodlim::neighborhood,
                         bandwidth = bandwidth)
      bandwidth <- sapply(nbh.list, function(nbh) nbh$bandwidth)
      tabS <- c(0, cumsum(tabulate(Sfactor))[-NS])
      neighbors <- unlist(lapply(1:NS, function(l) {
        nbh.list[[l]]$neighbors + tabS[l]
      }), use.names = FALSE)
    }
    response <- response[neighbors, , drop = FALSE]
    # if (weighted == TRUE)
    #   weights <- weights[neighbors]
  }

  switch(cotype, {
    size.strata <- NROW(response)
    NU <- 1
    N <- length(unique(response[,  "time"]))

  }, {
    size.strata <- tabulate(Sfactor)
    N <- NROW(response)
    NU <- length(size.strata)

  }, {
    size.strata <- nbh$size.nbh
    N <- sum(size.strata)
    NU <- nbh$nu

  }, {
    size.strata <- unlist(lapply(nbh.list, function(nbh) nbh$size.nbh),
                          use.names = FALSE)
    N <- sum(size.strata)

    n.unique.strata <- unlist(lapply(nbh.list, function(nbh) nbh$nu),
                              use.names = FALSE)
    NU <- sum(n.unique.strata)
  })

  if(min(size.strata)<xlimits){
    stop(paste("The minimum size of strata shall greater than ", xlimits, sep = ""))
  }


  X <- switch(cotype, {
    NULL
  }, {
    X <- data.frame(unique(strata[sorted, , drop = FALSE]))
    rownames(X) <- 1:NROW(X)
    X
  }, {
    X <- unlist(lapply(nbh.list, function(x) x$values), use.names = FALSE)
    X <- data.frame(X)
    colnames(X) <- colnames(NN)
    rownames(X) <- 1:NROW(X)
    X
  }, {
    D <- data.frame(unique(strata[sorted, , drop = FALSE]))
    D <- data.frame(D[rep(1:NS, n.unique.strata), , drop = FALSE])
    C <- data.frame(unlist(lapply(nbh.list, function(x) x$values),
                           use.names = FALSE))
    X <- cbind(D, C)
    colnames(X) <- c(colnames(strata), colnames(NN))
    rownames(X) <- 1:NROW(X)
    X
  }, {
    X = data.frame(pseudo = "pseudo")
    rownames(X) <- 1:NROW(X)
    X
  })

  continuous.predictors <- colnames(NN)
  discrete.predictors <- colnames(strata)
  model.matrix <- switch(cotype, {
    NULL
  }, strata, NN, cbind(strata, NN))[event.time.order, ,
                                    drop = FALSE]
  # event.history <- event.history[event.time.order, , drop = FALSE] ## sorted by time/id


  out <- list(time = response[, "time"], event = response[, "status"],#weights = weights,
              nn=N,nstrata = NU, stratas = X, size.strata=size.strata, # sorted by treatment
              covariate.type = cotype, xlevels = xlevels,
              discrete.predictors = discrete.predictors,
              continuous.predictors = continuous.predictors,
              model.response = event.history, model.matrix = model.matrix,
              # originalDataOrder = order(event.time.order),
              na.action = na.action,formula=formula)

  # event.history[ originalDataOrder,] is same with the original data
  if (cotype %in% c(3, 4)) out <- c(out, list(bandwidth = bandwidth))

  return(out)
}

#' @aliases prepformula
#' @noRd
#'
prepformula<- function (form1,form2) {

  # Call <- match.call()

  terms.labels <- union(attr(terms(form1),"term.labels"),attr(terms(form2),"term.labels"))
  if(length(terms.labels)>0){
    updated.terms <- paste(terms.labels,collapse=" + ")
    form1 <- as.formula(update(form1,paste(".~",updated.terms,collapse="")))

    form2 <- as.formula(update(form2,paste(".~",updated.terms,collapse="")))
  }
  return(list(form1=form1,form2=form2))

}


# prepData<- function (formula, data, weights, subset, na.action) {
#
#   Call <- match.call()
#   if (missing(formula))
#     stop("a formula argument is required")
#
#
#
#   indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(Call), nomatch = 0)
#   if (indx[1] == 0)
#     stop("A formula argument is required")
#   tform <- Call[c(1, indx)]
#   tform[[1L]] <- quote(stats::model.frame)
#
#
#   if (missing(data)) {tform$formula <- terms(formula, special)
#   }  else {tform$formula <-terms(formula,  data = data)}
#
#   mf <- eval(tform, parent.frame())
#   Terms <- terms(mf)
#   n <- nrow(mf)
#   Y <- model.response(mf)
#   if (!is.Surv(Y))
#     stop("Response must be a survival object")
#   id <- model.extract(mf, "id")
#
#   if (n == 0)
#     stop("No (non-missing) observations")
#   type <- attr(Y, "type")
#   if (type != "right" ) {
#     stop(paste("Cox model doesn't support \"", type, "\" survival data",
#                sep = ""))
#   }
#
#
#
#   if (length(attr(Terms, "variables")) > 2) {
#     ytemp <- innerterms(formula[1:2])
#     suppressWarnings(z <- as.numeric(ytemp))
#     ytemp <- ytemp[is.na(z)]
#     xtemp <- innerterms(formula[-2])
#     if (any(!is.na(match(xtemp, ytemp))))
#       warning("a variable appears on both the left and right sides of the formula")
#   }
#
#   dropterms <- NULL
#
#
#   xlevels <- .getXlevels(Terms, mf)
#   weights <- model.weights(mf)
#   has.id <- !(missing(id) || length(id) == 0)
#   has.rwt <- (!is.null(weights) && any(weights != floor(weights)))
#
#   if (has.id)
#     id <- as.factor(id)
#
#
#
#   contrast.arg <- NULL
#   attr(Terms, "intercept") <- 1
#
#   if (length(dropterms)) {
#     Terms2 <- Terms[-dropterms]
#     X <- model.matrix(Terms2, mf, constrasts.arg = contrast.arg)
#     temp <- attr(X, "assign")
#     shift <- sort(dropterms)
#     for (i in seq(along.with = shift)) temp <- temp + 1 *
#       (shift[i] <= temp)
#     attr(X, "assign") <- temp
#   }
#   else X <- model.matrix(Terms, mf, contrasts.arg = contrast.arg)
#   Xatt <- attributes(X)
#
#   adrop <- 0
#   xdrop <- Xatt$assign %in% adrop
#   X <- X[, !xdrop, drop = FALSE]
#   attr(X, "assign") <- Xatt$assign[!xdrop]
#   attr(X, "contrasts") <- Xatt$contrasts
#   offset <- model.offset(mf)
#   if (is.null(offset) || all(offset == 0)) {
#     offset <- rep(0, nrow(mf))
#     meanoffset <- 0
#   }
#   else if (any(!is.finite(offset) | !is.finite(exp(offset))))
#     stop("offsets must lead to a finite risk score")
#   else {
#     meanoffset <- mean(offset)
#     offset <- offset - meanoffset
#   }
#
#   if (!is.null(weights) && any(!is.finite(weights)))
#     stop("weights must be finite")
#
#
#   if (!all(is.finite(X)))
#     stop("data contains an infinite predictor")
#
#
#   return(list(time = Y[,1],event = Y[,2],X=X,weights=weights,offset = offset))
#
#
# }
