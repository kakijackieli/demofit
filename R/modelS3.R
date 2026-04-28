#' Common S3 methods for stochastic mortality models
#'
#' @description
#' Common S3 methods for objects returned by \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, \code{STARS()}, \code{ENS()}, \code{ENI()}, \code{CFMS()}, \code{CFM2S()}, and \code{CAES()}:
#'
#' coef(object) gives estimated parameter values.
#'
#' forecast::forecast(object) gives mortality forecasts smoothed by selected mortality curve (or ensemble interval forecasts for \code{ENI()}).
#'
#' plot(object) plots estimated parameter values, standardised residuals heatmap, data, and mortality forecasts smoothed by selected mortality curve (or ensemble interval forecasts for \code{ENI()}).
#'
#' residuals(object) gives standardised residuals on \eqn{ln(m_{x,t})}.
#'
#' simulate(object,nsim,seed) gives simulated future mortality rates (unsmoothed).
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
#' @method forecast ENI
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
#' @method plot ENI
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
#' @method simulate LCS
#' @method simulate RHS
#' @method simulate APCS
#' @method simulate CBDS
#' @method simulate CBDCS
#' @method simulate CBDQCS
#' @method simulate STARS
#'
NULL
