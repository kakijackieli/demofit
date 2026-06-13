#' Common S3 methods for stochastic mortality models
#'
#' @description
#' Common S3 methods for objects returned by \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, \code{STARS()}, \code{ENS()}, \code{ENI()}, \code{CFMS()}, \code{CFM2S()}, \code{CAES()}, \code{PLCS()}, \code{PRHS()}, \code{PAPCS()}, \code{PCBDS()}, \code{PCBDCS()}, \code{PCBDQCS()}, \code{PCFMS()}, \code{PCFM2S()}, \code{PCAES()}, and \code{LCLS()}.
#'
#' coef(object) gives estimated parameter values.
#'
#' forecast::forecast(object) gives mortality forecasts (which = 1 for smoothed (default); which = 2 for raw) (or ensemble interval forecasts for \code{ENI()}).
#'
#' plot(object) plots estimated parameter values, standardised residuals heatmap, data, and mortality forecasts smoothed by selected mortality curve (or ensemble interval forecasts for \code{ENI()}).
#'
#' residuals(object) gives standardised residuals on \eqn{ln(m_{x,t})} or \eqn{D_{x,t}}.
#'
#' simulate(object,nsim,seed) gives simulated future mortality rates (raw).
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
#' @method coef PLCS
#' @method coef PRHS
#' @method coef PAPCS
#' @method coef PCBDS
#' @method coef PCBDCS
#' @method coef PCBDQCS
#' @method coef PCFMS
#' @method coef PCFM2S
#' @method coef PCAES
#' @method coef LCLS
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
#' @method forecast PLCS
#' @method forecast PRHS
#' @method forecast PAPCS
#' @method forecast PCBDS
#' @method forecast PCBDCS
#' @method forecast PCBDQCS
#' @method forecast PCFMS
#' @method forecast PCFM2S
#' @method forecast PCAES
#' @method forecast LCLS
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
#' @method plot PLCS
#' @method plot PRHS
#' @method plot PAPCS
#' @method plot PCBDS
#' @method plot PCBDCS
#' @method plot PCBDQCS
#' @method plot PCFMS
#' @method plot PCFM2S
#' @method plot PCAES
#' @method plot LCLS
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
#' @method residuals PLCS
#' @method residuals PRHS
#' @method residuals PAPCS
#' @method residuals PCBDS
#' @method residuals PCBDCS
#' @method residuals PCBDQCS
#' @method residuals PCFMS
#' @method residuals PCFM2S
#' @method residuals PCAES
#' @method residuals LCLS
#'
#' @method simulate LCS
#' @method simulate RHS
#' @method simulate APCS
#' @method simulate CBDS
#' @method simulate CBDCS
#' @method simulate CBDQCS
#' @method simulate STARS
#' @method simulate LCLS
#'
NULL
