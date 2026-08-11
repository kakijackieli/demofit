#' Poisson mortality model fitting
#'
#' Fits and forecasts mortality rates using different stochastic mortality models with Poisson assumption.
#'
#' @param x vector of ages.
#' @param D matrix of death counts (rows as years and columns as ages).
#' @param E matrix of mid-year exposures (rows as years and columns as ages).
#' @param model name of stochastic mortality model (including LC, RH, APC, CBD, CBDC, and CBDQC).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive) (optional; if missing, no smoothing is performed).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1).
#'
#' @details
#' See \code{PLCS()}, \code{PRHS()}, \code{PAPCS()}, \code{PCBDS()}, \code{PCBDCS()}, and \code{PCBDQCS()} for more details of different stochastic mortality models with Poisson assumption.
#'
#' @return 
#' An object of class based on selected stochastic mortality model with associated S3 methods coef, forecast, plot, and residuals.
#'
#' @examples
#' x <- 60:89
#' a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962,
#' -3.8885,-3.7896,-3.6853,-3.5737,-3.4728,-3.3718,-3.2586,-3.1474,-3.0371,-2.9206,
#' -2.7998,-2.6845,-2.5653,-2.4581,-2.3367,-2.2159,-2.1017,-1.9941,-1.8821, -1.7697)
#' b <- c(0.0283,0.0321,0.0335,0.0336,0.0341,0.0358,0.0368,0.0403,0.0392,0.0395,
#' 0.0396,0.0399,0.0397,0.0386,0.039,0.0375,0.0367,0.0368,0.035,0.0354,
#' 0.0336,0.0323,0.0313,0.0295,0.0282,0.0265,0.024,0.0226,0.0219,0.0183)
#' k <- c(12.11,10.69,11.18,9.64,9.35,8.21,6.89,5.74,4.56,3.6,
#' 3.27,2.04,1.11,-0.44,-1.05,-1.03,-1.84,-2.9,-4.03,-4.12,
#' -5.18,-5.64,-6,-6.51,-6.91,-6.9,-8.32,-8.53,-9.69,-9.31)
#' set.seed(123)
#' M <- exp(outer(k,b)+matrix(a,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.035))
#' E <- matrix(c(107788,108036,107481,106552,104608,100104,95803,91345,84980,79557,
#' 75146,70559,65972,60898,55623,50522,47430,45895,41443,34774,
#' 30531,27754,25105,22271,19437,16888,14458,12146,10038,7994),30,30,byrow=TRUE)
#' D <- round(E*M)
#' fit <- PFCS(x=x,D=D,E=E,model="LC",curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
PFCS <- function(x,D,E,model=c("LC","RH","APC","CBD","CBDC","CBDQC"),curve=NULL,h=10,jumpoff=1) {
model <- tryCatch(match.arg(model),error = function(e) { stop("invalid model choice") })
if (!is.null(curve)) { curve <- tryCatch(match.arg(curve,choices=.curves),error = function(e) { stop("invalid curve choice") }) }
tryCatch({
fit <- switch(model,
LC = PLCS(x,D,E,curve,h,jumpoff),
RH = PRHS(x,D,E,curve,h,jumpoff),
APC = PAPCS(x,D,E,curve,h,jumpoff),
CBD = PCBDS(x,D,E,curve,h,jumpoff),
CBDC = PCBDCS(x,D,E,curve,h,jumpoff),
CBDQC = PCBDQCS(x,D,E,curve,h,jumpoff)
)
invisible(fit)
}, error = function(e) { stop(paste0("model fitting and forecasting is unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}
