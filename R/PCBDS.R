#' Poisson CBD model
#'
#' Fits and forecasts mortality rates using CBD model with Poisson assumption.
#'
#' @param x vector of ages.
#' @param D matrix of death counts (rows as years and columns as ages).
#' @param E matrix of mid-year exposures (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, perks, weibull, beard, martinelle, thatcher, gompertz2, makeham2, perks2, weibull2, beard2, martinelle2, thatcher2, where first 7 curves' parameters are unconstrained and last 7 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The CBD (M5) model with Poisson assumption is specified as 
#' 
#' \eqn{ln(m_{x,t}) = \kappa_{1,t} + \kappa_{2,t} (x-\bar{x})} and \eqn{D_{x,t} ~ Poisson(E_{x,t} m_{x,t})}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_{1,t}} and \eqn{\kappa_{2,t}}. It is designed for ages 50-90.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted lm sd simulate rnorm
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class PCBDS with associated S3 methods coef, forecast (which = 1 for smoothed (default); which = 2 for raw), plot, and residuals.
#'
#' @references
#' Cairns, A.J.G., Blake, D., and Dowd, K. (2006). A two-factor model for stochastic mortality with parameter uncertainty: Theory and calibration. Journal of Risk and Insurance, 73(4), 687-718.
#'
#' @examples
#' x <- 60:89
#' k1 <- -2.97-0.0245*(0:29)
#' k2 <- 0.101+0.000345*(0:29)
#' set.seed(123)
#' M <- exp(matrix(k1,nrow=30,ncol=30,byrow=FALSE)+outer(k2,(x-mean(x)))+rnorm(900,0,0.035))
#' E <- matrix(c(107788,108036,107481,106552,104608,100104,95803,91345,84980,79557,
#' 75146,70559,65972,60898,55623,50522,47430,45895,41443,34774,
#' 30531,27754,25105,22271,19437,16888,14458,12146,10038,7994),30,30,byrow=TRUE)
#' D <- round(E*M)
#' fit <- PCBDS(x=x,D=D,E=E,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
PCBDS <- function(x,D,E,curve=c("gompertz","makeham","perks","weibull","beard","martinelle","thatcher","gompertz2","makeham2","perks2","weibull2","beard2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(D)||!is.numeric(E)) { stop("x and D and E must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(D)||!is.matrix(E)||nrow(D)!=nrow(E)) stop("D and E must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(D)||length(x)!=ncol(E)) stop("the number of ages must match the number of columns of D and E")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if ((min(x)<50)||(max(x)>90)) { stop("x must be in between 50 and 90") }
if (any(D<=0)||any(E<=0)) { stop("all D and E values must be positive") }
if (nrow(D)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
M <- D/E
nr <- nrow(M); nc <- ncol(M)
k1 <- numeric(); k2 <- numeric()
xx <- x-mean(x); for (i in 1:nr) { yy <- log(M[i,]); k1[i] <- as.numeric(lm(yy~xx)$coef[1]); k2[i] <- as.numeric(lm(yy~xx)$coef[2]) }
olde <- 1000000000; tol <- 1e-8
for (z in 1:200) {
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); xm <- matrix(x-mean(x),nr,nc,byrow=TRUE); Emat <- E*exp(k1mat+k2mat*xm)
k1 <- k1+rowSums(D-Emat)/rowSums(Emat)
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); xm <- matrix(x-mean(x),nr,nc,byrow=TRUE); Emat <- E*exp(k1mat+k2mat*xm)
k2 <- k2+rowSums((D-Emat)*xm)/rowSums(Emat*xm^2)
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); xm <- matrix(x-mean(x),nr,nc,byrow=TRUE); Emat <- E*exp(k1mat+k2mat*xm)
newe <- sum(D*log(D/Emat)-D+Emat)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
dis <- sum(2*(D*log(D/Emat)-D+Emat))/(nr*nc-length(k1)-length(k2))
res <- sign(D-Emat)*sqrt(2*(D*log(D/Emat)-D+Emat)/dis)
k1f <- suppressMessages(forecast(auto.arima(tsclean(k1)),h=h)$mean)
k2f <- suppressMessages(forecast(auto.arima(tsclean(k2)),h=h)$mean)
Mf <- array(NA,c(h,nc)); Mfs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- exp(k1f[i]+k2f[i]*(x[j]-mean(x))) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- M[nr,j]*exp(k1f[i]-k1[nr]+(k2f[i]-k2[nr])*(x[j]-mean(x))) }}}
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,D=D,E=E,M=M,h=h,jumpoff=jumpoff,kappa1=k1,kappa2=k2,standardresiduals=res,dispersion=dis,forecast=Mf,smoothforecast=Mfs),
class="PCBDS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.PCBDS <- function(object,...) {
list(kappa1=object$kappa1,kappa2=object$kappa2)
}

#' @export
forecast.PCBDS <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { object$smoothforecast }
else if (which==2) { object$forecast }
}

#' @export
plot.PCBDS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(2,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(c(1:nrow(x$M)),x$kappa1,xlab="year",ylab="kappa1",bty="n",type="l")
plot(c(1:nrow(x$M)),x$kappa2,xlab="year",ylab="kappa2",bty="n",type="l")
colband <- colorRampPalette(c("black","grey","white"))
image(x=x$x,y=c(1:nrow(x$M)),z=t(x$standardresiduals),col=colband(6),xlab="age",ylab="year",main="standardised residuals",cex.main=1,font.main=1)
plot(x$x,log(x$M[1,]),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M),log(x$smoothforecast)),max(log(x$M),log(x$smoothforecast))))
points(x$x,log(x$M[nrow(x$M),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
residuals.PCBDS <- function(object,...) {
object$standardresiduals
}
