prepData<- function (formula, data, weights, subset, na.action) {
  
  Call <- match.call()
  if (missing(formula)) 
    stop("a formula argument is required")
  
  
  
  indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  tform <- Call[c(1, indx)]
  tform[[1L]] <- quote(stats::model.frame)
  
  
  if (missing(data)) {tform$formula <- terms(formula, special)
  }  else {tform$formula <-terms(formula,  data = data)}
  
  mf <- eval(tform, parent.frame())
  Terms <- terms(mf)
  n <- nrow(mf)
  Y <- model.response(mf)
  if (!is.Surv(Y)) 
    stop("Response must be a survival object")
  id <- model.extract(mf, "id")
  
  if (n == 0) 
    stop("No (non-missing) observations")
  type <- attr(Y, "type")
  if (type != "right" ) {
    stop(paste("Cox model doesn't support \"", type, "\" survival data", 
               sep = ""))
  }
  
  
  
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- innerterms(formula[1:2])
    suppressWarnings(z <- as.numeric(ytemp))
    ytemp <- ytemp[is.na(z)]
    xtemp <- innerterms(formula[-2])
    if (any(!is.na(match(xtemp, ytemp)))) 
      warning("a variable appears on both the left and right sides of the formula")
  }
  
  dropterms <- NULL
  
  
  xlevels <- .getXlevels(Terms, mf)
  weights <- model.weights(mf)
  has.id <- !(missing(id) || length(id) == 0)
  has.rwt <- (!is.null(weights) && any(weights != floor(weights)))
  
  if (has.id) 
    id <- as.factor(id)
  
  
  
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  
  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts.arg = contrast.arg)
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    for (i in seq(along.with = shift)) temp <- temp + 1 * 
      (shift[i] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts.arg = contrast.arg)
  Xatt <- attributes(X)
  
  adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) || all(offset == 0)) {
    offset <- rep(0, nrow(mf))
    meanoffset <- 0
  }
  else if (any(!is.finite(offset) | !is.finite(exp(offset)))) 
    stop("offsets must lead to a finite risk score")
  else {
    meanoffset <- mean(offset)
    offset <- offset - meanoffset
  }
  
  if (!is.null(weights) && any(!is.finite(weights))) 
    stop("weights must be finite")
  
  
  if (!all(is.finite(X))) 
    stop("data contains an infinite predictor")
  
  
  return(list(time = Y[,1],event = Y[,2],X=X,weights=weights,offset = offset))
  
  
}