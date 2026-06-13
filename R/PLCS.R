#' Poisson Lee-Carter model
#'
#' Fits and forecasts mortality rates using Lee-Carter model with Poisson assumption.
#'
#' @param x vector of ages.
#' @param D matrix of death counts (rows as years and columns as ages).
#' @param E matrix of mid-year exposures (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The Lee-Carter (LC) model with Poisson assumption is specified as 
#' 
#' \eqn{ln(m_{x,t}) = \alpha_x + \beta_x \kappa_t} and \eqn{D_{x,t} ~ Poisson(E_{x,t} m_{x,t})}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_t}. Constraints include sum of \eqn{\beta_x} is one and sum of \eqn{\kappa_t} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd simulate rnorm
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class PLCS with associated S3 methods coef, forecast (which = 1 for smoothed (default); which = 2 for raw), plot, and residuals.
#'
#' @references
#' Renshaw, A.E. and Haberman, S. (2006). A cohort-based extension to the Lee–Carter model for mortality reduction factors. Insurance: Mathematics and Economics, 38(3), 556-570.
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
#' fit <- PLCS(x=x,D=D,E=E,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
PLCS <- function(x,D,E,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(D)||!is.numeric(E)) { stop("x and D and E must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(D)||!is.matrix(E)||nrow(D)!=nrow(E)) stop("D and E must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(D)||length(x)!=ncol(E)) stop("the number of ages must match the number of columns of D and E")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(D<=0)||any(E<=0)) { stop("all D and E values must be positive") }
if (nrow(D)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
M <- D/E
nr <- nrow(M); nc <- ncol(M)
PCA <- prcomp(log(M),center=TRUE,scale=FALSE)
a <- numeric(); for (j in 1:nc) { a[j] <- mean(log(M[,j])) }
b <- PCA$rotation[,1]/sum(PCA$rotation[,1])
k <- PCA$x[,1]*sum(PCA$rotation[,1])
olde <- 1000000000; tol <- 1e-8
for (z in 1:200) {
amat <- matrix(a,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); Emat <- E*exp(amat+bmat*kmat)
a <- a+colSums(D-Emat)/colSums(Emat)
amat <- matrix(a,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); Emat <- E*exp(amat+bmat*kmat)
b <- b+colSums((D-Emat)*kmat)/colSums(Emat*kmat^2); b <- b/sum(b)
amat <- matrix(a,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); Emat <- E*exp(amat+bmat*kmat)
k <- k+rowSums((D-Emat)*bmat)/rowSums(Emat*bmat^2); k <- k-mean(k)
amat <- matrix(a,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); Emat <- E*exp(amat+bmat*kmat)
newe <- sum(D*log(D/Emat)-D+Emat)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
dis <- sum(2*(D*log(D/Emat)-D+Emat))/(nr*nc-length(a)-length(b)-length(k)+2)
res <- sign(D-Emat)*sqrt(2*(D*log(D/Emat)-D+Emat)/dis)
kf <- suppressMessages(forecast(auto.arima(tsclean(k)),h=h)$mean)
Mf <- array(NA,c(h,nc)); Mfs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- exp(a[j]+b[j]*kf[i]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- M[nr,j]*exp(b[j]*(kf[i]-k[nr])) }}}
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,D=D,E=E,M=M,h=h,jumpoff=jumpoff,alpha=a,beta=b,kappa=k,standardresiduals=res,dispersion=dis,forecast=Mf,smoothforecast=Mfs),
class="PLCS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.PLCS <- function(object,...) {
list(alpha=object$alpha,beta=object$beta,kappa=object$kappa)
}

#' @export
forecast.PLCS <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { object$smoothforecast }
else if (which==2) { object$forecast }
}

#' @export
plot.PLCS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(3,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha,xlab="age",ylab="alpha",bty="n",type="l")
plot(x$x,x$beta,xlab="age",ylab="beta",bty="n",type="l")
plot(c(1:nrow(x$M)),x$kappa,xlab="year",ylab="kappa",bty="n",type="l")
colband <- colorRampPalette(c("black","grey","white"))
image(x=x$x,y=c(1:nrow(x$M)),z=t(x$standardresiduals),col=colband(6),xlab="age",ylab="year",main="standardised residuals",cex.main=1,font.main=1)
plot(x$x,log(x$M[1,]),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M),log(x$smoothforecast)),max(log(x$M),log(x$smoothforecast))))
points(x$x,log(x$M[nrow(x$M),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
residuals.PLCS <- function(object,...) {
object$standardresiduals
}
