#################################################################
## Accessor functions for CopulaSCR objects
#################################################################


#' Extract Predicted Values
#'
#' Extracts the principal predicted values from an object returned by a
#' prediction function in the CopulaSCR package. This function provides a
#' stable interface for accessing prediction results without requiring users
#' to know the internal list structure of the prediction object.
#'
#' @param object A prediction object returned by \code{predjsurv()},
#'   \code{predcsurv()}, \code{predchr()}, \code{predict.scrsurv()}, or
#'   \code{predict.mscr()}.
#' @param event For marginal survival predictions, a character string
#'   specifying the event type to extract. Possible values are
#'   \code{"both"}, \code{"nonterminal"}, and \code{"terminal"}.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return The principal predicted values. The structure of the returned
#'   object depends on the prediction type.
#'
#' @seealso \code{\link{predictscr}},
#'   \code{\link{predict.scrsurv}},
#'   \code{\link{predict.mscr}}
#'
#' @export
predicted_values <- function(object, ...) {
  UseMethod("predicted_values")
}


#' @rdname predicted_values
#' @export
predicted_values.predjsurv <- function(object, ...) {
  object$jsurv
}


#' @rdname predicted_values
#' @export
predicted_values.predcsurv <- function(object, ...) {
  object$csurv
}


#' @rdname predicted_values
#' @export
predicted_values.predchr <- function(object, ...) {
  object$chr
}


#' @rdname predicted_values
#' @export
predicted_values.predmsurv <- function(
    object,
    event = c("both", "nonterminal", "terminal"),
    ...) {
  
  event <- match.arg(event)
  
  if (event == "nonterminal") {
    if (is.null(object$msurv1)) {
      stop(
        "Nonterminal-event survival predictions are not available.",
        call. = FALSE
      )
    }
    
    return(object$msurv1)
  }
  
  if (event == "terminal") {
    if (is.null(object$msurv2)) {
      stop(
        "Terminal-event survival predictions are not available.",
        call. = FALSE
      )
    }
    
    return(object$msurv2)
  }
  
  has.nonterminal <- !is.null(object$msurv1)
  has.terminal <- !is.null(object$msurv2)
  
  if (has.nonterminal && has.terminal) {
    return(
      list(
        nonterminal = object$msurv1,
        terminal = object$msurv2
      )
    )
  }
  
  if (has.nonterminal) {
    return(object$msurv1)
  }
  
  if (has.terminal) {
    return(object$msurv2)
  }
  
  stop(
    "No marginal survival predictions are available.",
    call. = FALSE
  )
}


#' @rdname predicted_values
#' @export
predicted_values.predmscr <- function(object, ...) {
  
  type <- attr(object, "type", exact = TRUE)
  
  if (is.null(type)) {
    stop(
      "The prediction type is not stored in the object.",
      call. = FALSE
    )
  }
  
  if (identical(type, "msurv")) {
    return(
      lapply(
        unclass(object),
        predicted_values
      )
    )
  }
  
  if (!type %in% names(object)) {
    stop(
      "The requested prediction component is not available.",
      call. = FALSE
    )
  }
  
  value <- object[[type]]
  
  # For a single subject, the principal prediction is already stored
  # directly in this component.
  if (!is.list(value)) {
    return(value)
  }
  
  # For multiple subjects, each element may contain a subject-specific
  # prediction record.
  is.prediction.record <- vapply(
    value,
    function(x) {
      is.list(x) && type %in% names(x)
    },
    logical(1)
  )
  
  if (length(value) > 0L &&
      all(is.prediction.record)) {
    
    value <- lapply(
      value,
      function(x) x[[type]]
    )
    
    if (all(lengths(value) == 1L)) {
      return(
        unlist(value, use.names = FALSE)
      )
    }
  }
  
  value
}



#' Extract Association Estimates
#'
#' Extracts the estimated Kendall's tau values or copula parameters and their
#' standard errors from an object returned by \code{scrassonp()}. This function
#' avoids direct access to the internal components of a fitted
#' \code{"scrassonp"} object.
#'
#' @param object An object of class \code{"scrassonp"}.
#' @param scale A character string specifying the association measure to
#'   extract. Possible values are \code{"tau"} (the default) and
#'   \code{"copula"}.
#' @param ... Additional arguments; currently unused.
#'
#' @return A data frame containing the estimates and their standard errors.
#'   Standard errors are reported as \code{NA} when they were not computed
#'   during model fitting.
#'
#' @seealso \code{\link{scrassonp}},
#'   \code{\link{summary.scrassonp}}
#'
#' @export
association_estimates <- function(object, ...) {
  UseMethod("association_estimates")
}


#' @rdname association_estimates
#' @export
association_estimates.scrassonp <- function(
    object,
    scale = c("tau", "copula"),
    ...) {
  
  scale <- match.arg(scale)
  
  if (scale == "tau") {
    estimate <- as.numeric(object$tau)
    std.error <- object$tau.se
  } else {
    estimate <- as.numeric(object$copulaparam)
    std.error <- object$param.se
  }
  
  if (is.null(std.error) ||
      length(std.error) != length(estimate)) {
    std.error <- rep(NA_real_, length(estimate))
  } else {
    std.error <- as.numeric(std.error)
  }
  
  out <- data.frame(
    estimate = estimate,
    std.error = std.error
  )
  
  if (length(estimate) == 1L) {
    rownames(out) <- "overall"
  } else {
    rownames(out) <- paste0(
      "strata",
      seq_along(estimate)
    )
  }
  
  out
}



#' Extract a Marginal Fit
#'
#' Extracts the fitted marginal model for one intermediate event from an
#' object returned by \code{mscr()}. This function avoids direct access to the
#' internal \code{mar.fits} component of an \code{"mscr"} object.
#'
#' @param object An object of class \code{"mscr"}.
#' @param event A single integer or character value identifying the
#'   intermediate event. An integer selects the event by position, whereas
#'   a character value selects it by name.
#' @param ... Additional arguments; currently unused.
#'
#' @return The corresponding fitted object of class \code{"scrsurv"}.
#'
#' @seealso \code{\link{mscr}}, \code{\link{scrsurv}}
#'
#' @export
marginal_fit <- function(object, ...) {
  UseMethod("marginal_fit")
}


#' @rdname marginal_fit
#' @export
marginal_fit.mscr <- function(object, event, ...) {
  
  if (missing(event) || length(event) != 1L) {
    stop(
      "'event' must identify exactly one intermediate event.",
      call. = FALSE
    )
  }
  
  fits <- object$mar.fits
  
  if (is.character(event)) {
    if (is.null(names(fits)) ||
        !event %in% names(fits)) {
      stop(
        "The specified intermediate-event name was not found.",
        call. = FALSE
      )
    }
    
    return(fits[[event]])
  }
  
  if (!is.numeric(event) ||
      is.na(event) ||
      event != as.integer(event) ||
      event < 1L ||
      event > length(fits)) {
    stop(
      "'event' must be a valid intermediate-event index or name.",
      call. = FALSE
    )
  }
  
  fits[[as.integer(event)]]
}



#' Extract the Terminal Survival Function
#'
#' Extracts the estimated marginal survival function of the terminal event
#' from an object returned by \code{mscr()}. This function avoids direct access
#' to the internal \code{S_D} component of the fitted object.
#'
#' @param object An object of class \code{"mscr"}.
#' @param ... Additional arguments; currently unused.
#'
#' @return The estimated terminal-event survival function.
#'
#' @seealso \code{\link{mscr}}, \code{\link{dyBS}}
#'
#' @export
terminal_survival <- function(object, ...) {
  UseMethod("terminal_survival")
}


#' @rdname terminal_survival
#' @export
terminal_survival.mscr <- function(object, ...) {
  object$S_D
}



#' Extract the Integrated Brier Score
#'
#' Extracts the integrated Brier score or the integrated relative Brier score
#' from an object returned by \code{dyBS()}.
#'
#' @param object An object of class \code{"dyBS"}.
#' @param relative Logical; whether to return the integrated relative Brier
#'   score. This quantity is available only when \code{dyBS()} was called with
#'   \code{reference = TRUE}.
#' @param ... Additional arguments; currently unused.
#'
#' @return A numeric value containing the requested integrated score.
#'
#' @seealso \code{\link{dyBS}}
#'
#' @export
integrated_brier_score <- function(object, ...) {
  UseMethod("integrated_brier_score")
}


#' @rdname integrated_brier_score
#' @export
integrated_brier_score.dyBS <- function(
    object,
    relative = FALSE,
    ...) {
  
  if (isTRUE(relative)) {
    if (is.null(object$IRBS)) {
      stop(
        "The integrated relative Brier score is not available because ",
        "the object was created with 'reference = FALSE'.",
        call. = FALSE
      )
    }
    
    return(object$IRBS)
  }
  
  object$IBS
}
#' Extract Brier Scores
#'
#' Extracts the Brier-score curve from an object returned by \code{dyBS()}.
#'
#' @param object An object of class \code{"dyBS"}.
#' @param ... Additional arguments; currently unused.
#'
#' @return A data frame containing the evaluation times and estimated
#'   Brier scores.
#'
#' @seealso \code{\link{dyBS}}, \code{\link{integrated_brier_score}}
#'
#' @export
brier_scores <- function(object, ...) {
  UseMethod("brier_scores")
}

#' @rdname brier_scores
#' @export
brier_scores.dyBS <- function(object, ...) {
  object$BS
}

#' Extract Pseudo-Log-Likelihood from an 'mscr' Object
#'
#' Extracts the pseudo-log-likelihood value from a fitted
#' \code{"mscr"} object.
#'
#' @param object An object of class \code{"mscr"}.
#' @param ... Further arguments; currently unused.
#'
#' @return The pseudo-log-likelihood value.
#'
#' @seealso \code{\link{mscr}}
#'
#' @export
logLik.mscr <- function(object, ...) {
  object$logLik.alpha
}