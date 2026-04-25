#' Common S3 methods for stochastic mortality models
#'
#' @description
#' Common S3 methods for objects returned by \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, and \code{STARS()}:
#'
#' coef(object) gives estimated parameter values.
#'
#' forecast::forecast(object) gives mortality forecasts smoothed by selected mortality curve.
#'
#' plot(object) plots estimated parameter values, standardised residuals heatmap, data, and mortality forecasts smoothed by selected mortality curve.
#'
#' residuals(object) gives standardised residuals on \eqn{ln(m_{x,t})}.
#'
#' @name modelS3
#'
#' @method coef LCS
#' @method coef RHS
#' @method coef APCS
#' @method coef CBDS
#' @method coef CBDCS
#' @method coef CBDQCS
#' @method coef STARS
#' @method coef CFMS
#' @method coef CFM2S
#' @method coef CAES
#'
#' @method forecast LCS
#' @method forecast RHS
#' @method forecast APCS
#' @method forecast CBDS
#' @method forecast CBDCS
#' @method forecast CBDQCS
#' @method forecast STARS
#' @method forecast ENS
#' @method forecast CFMS
#' @method forecast CFM2S
#' @method forecast CAES
#'
#' @method plot LCS
#' @method plot RHS
#' @method plot APCS
#' @method plot CBDS
#' @method plot CBDCS
#' @method plot CBDQCS
#' @method plot STARS
#' @method plot ENS
#' @method plot CFMS
#' @method plot CFM2S
#' @method plot CAES
#'
#' @method residuals LCS
#' @method residuals RHS
#' @method residuals APCS
#' @method residuals CBDS
#' @method residuals CBDCS
#' @method residuals CBDQCS
#' @method residuals STARS
#' @method residuals CFMS
#' @method residuals CFM2S
#' @method residuals CAES
#'
NULL
