#' Mortality ensemble
#'
#' Fits and forecasts mortality rates using mortality ensemble.
#'
#' @param x vector of ages.
#' @param M matrix of mortality rates (rows as years and columns as ages).
#' @param wm vector of weights for LC, RH, APC, M5, M6, M7, and STAR models (default = 1/7).
#' @param curve name of mortality curve for smoothing ensemble mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1).
#'
#' @details
#' Ensemble forecast is obtained as a weighted average of forecasts from individual stochastic mortality models. See \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, and \code{STARS()} for more details of different stochastic mortality models.
#'
#' @return 
#' An object of class ENS with associated S3 methods forecast and plot.
#'
#' @references
#' Li, J. (2023). A model stacking approach for forecasting mortality. North American Actuarial Journal, 27(3), 530-545.
#'
#' @examples
#' x <- 60:69
#' a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
#' b <- c(0.0801,0.0909,0.0948,0.0951,0.0965,0.1014,0.1042,0.1141,0.1110,0.1118)
#' k <- c(12.11,10.69,11.18,9.64,9.35,8.21,6.89,5.74,4.56,3.60,
#' 3.27,2.04,1.11,-0.44,-1.05,-1.03,-1.84,-2.90,-4.03,-4.12,
#' -5.18,-5.64,-6.00,-6.51,-6.91,-6.90,-8.32,-8.53,-9.69,-9.31)
#' set.seed(123)
#' M <- exp(outer(k,b)+matrix(a,nrow=30,ncol=10,byrow=TRUE)+rnorm(300,0,0.035))
#' fit <- ENS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' forecast::forecast(fit)
#' plot(fit)
#'
#' @export
ENS <- function(x,M,wm=rep(1/7,7),curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (length(wm)!=7) { stop("wm must have a length of 7") }
if (any(wm<0)) { stop("wm must be non-negative") }
if (abs(sum(wm)-1)>1e-8) { stop("wm must sum to 1") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
fit1 <- LCS(x,M,curve,h,jumpoff)
fit2 <- RHS(x,M,curve,h,jumpoff)
fit3 <- APCS(x,M,curve,h,jumpoff)
fit4 <- CBDS(x,M,curve,h,jumpoff)
fit5 <- CBDCS(x,M,curve,h,jumpoff)
fit6 <- CBDQCS(x,M,curve,h,jumpoff)
fit7 <- STARS(x,M,curve,h)
Mf <- wm[1]*fit1$forecast+wm[2]*fit2$forecast+wm[3]*fit3$forecast+wm[4]*fit4$forecast+wm[5]*fit5$forecast+wm[6]*fit6$forecast+wm[7]*fit7$forecast
Mfs <- array(NA,c(h,ncol(M)))
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M=M,wm=wm,h=h,jumpoff=jumpoff,forecast=Mf,smoothforecast=Mfs),
class="ENS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting is unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
forecast.ENS <- function(object,...) {
object$smoothforecast
}

#' @export
plot.ENS <- function(x,...) {
plot(x$x,log(x$M[1,]),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M),log(x$smoothforecast)),max(log(x$M),log(x$smoothforecast))))
points(x$x,log(x$M[nrow(x$M),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}
