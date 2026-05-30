#' Example semi-competing risks dataset with a single intermediate event
#'
#' A simulated dataset for semi-competing risks analysis with a single
#' intermediate event and one terminal event.
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{T1}{Observed time for the intermediate event.}
#'   \item{T2}{Observed time for the terminal event.}
#'   \item{event1}{Event indicator for the intermediate event.}
#'   \item{event2}{Event indicator for the terminal event.}
#' }
#' @usage data(SCRdata)
"SCRdata"

#' Example semi-competing risks dataset by treatment group
#'
#' A simulated dataset for semi-competing risks analysis with a single
#' intermediate event, one terminal event, and a treatment variable.
#'
#' @format A data frame with 500 rows and 5 variables:
#' \describe{
#'   \item{T1}{Observed time for the intermediate event.}
#'   \item{T2}{Observed time for the terminal event.}
#'   \item{event1}{Event indicator for the intermediate event.}
#'   \item{event2}{Event indicator for the terminal event.}
#'   \item{tr}{Treatment group indicator.}
#' }
#' @usage data(SCRdata_by_tr)
"SCRdata_by_tr"

#' Example semi-competing risks dataset with multiple intermediate events
#'
#' A simulated dataset for semi-competing risks analysis with three
#' intermediate events and one terminal event.
#'
#' @format An object of class \code{mscrData}, implemented as a list with
#' five components:
#' \describe{
#'   \item{T1}{A numeric matrix with 100 rows and 3 columns, containing observed times for three intermediate events.}
#'   \item{T2}{A numeric vector of length 100, containing observed terminal event times.}
#'   \item{event1}{A numeric matrix with 100 rows and 3 columns, containing event indicators for the three intermediate events.}
#'   \item{event2}{A numeric vector of length 100, containing terminal event indicators.}
#'   \item{call}{The function call used to generate the data.}
#' }
#' @usage data(mSCRdata)
"mSCRdata"