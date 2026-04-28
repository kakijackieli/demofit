#' Age-period-cohort model
#'
#' Fits and forecasts mortality rates using age-period-cohort model.
#'
#' @param x vector of ages.
#' @param M matrix of mortality rates (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The age-period-cohort (APC) model is specified as 
#' 
#' \eqn{ln(m_{x,t}) = \alpha_x + \kappa_t + \gamma_{t-x} + \epsilon_{x,t}}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_t} and \eqn{\gamma_c}. Constraints include sum of \eqn{\kappa_t} is zero and and sum of \eqn{\gamma_c} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd simulate rnorm
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class APCS with associated S3 methods coef, forecast, plot, residuals, and simulate (nsim for setting number of simulations; seed for initialising random number generator).
#'
#' @references
#' Bray, I. (2002). Application of Markov chain Monte Carlo methods to projecting cancer incidence and mortality. Journal of the Royal Statistical Society Series C, 51(2), 151-164.
#'
#' @examples
#' x <- 60:89
#' a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962,
#' -3.8885,-3.7896,-3.6853,-3.5737,-3.4728,-3.3718,-3.2586,-3.1474,-3.0371,-2.9206,
#' -2.7998,-2.6845,-2.5653,-2.4581,-2.3367,-2.2159,-2.1017,-1.9941,-1.8821, -1.7697)
#' b <- rep(1/30,30)
#' k <- c(12.11,10.69,11.18,9.64,9.35,8.21,6.89,5.74,4.56,3.6,
#' 3.27,2.04,1.11,-0.44,-1.05,-1.03,-1.84,-2.9,-4.03,-4.12,
#' -5.18,-5.64,-6,-6.51,-6.91,-6.9,-8.32,-8.53,-9.69,-9.31)
#' set.seed(123)
#' M <- exp(outer(k,b)+matrix(a,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.035))
#' fit <- APCS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
APCS <- function(x,M,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(M)) { stop("x and M must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M)) stop("M must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M)) stop("the number of ages must match the number of columns of M")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(M<=0)) { stop("all M values must be positive") }
if (nrow(M)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
nr <- nrow(M); nc <- ncol(M)
PCA <- prcomp(log(M),center=TRUE,scale=FALSE)
a <- numeric(); for (j in 1:nc) { a[j] <- mean(log(M[,j])) }
k <- PCA$x[,1]*sum(PCA$rotation[,1])
g <- rep(0,nr+nc-1)
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
amat <- matrix(a,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]
a <- a+colSums(log(M)-amat-kmat-gmat)/nr
amat <- matrix(a,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]
k <- k+rowSums(log(M)-amat-kmat-gmat)/nc; k <- k-mean(k)
amat <- matrix(a,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]
g <- g+as.numeric(tapply(as.vector(log(M)-amat-kmat-gmat),as.vector(row(M)-col(M)+nc),mean)); g <- g-mean(g)
amat <- matrix(a,nr,nc,byrow=TRUE); kmat <- matrix(k,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]
newe <- sum((log(M)-amat-kmat-gmat)^2)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
res <- array(NA,c(nr,nc))
for (i in 1:nr) { for (j in 1:nc) { res[i,j] <- log(M[i,j])-a[j]-k[i]-g[i-j+nc] }}
res <- (res-mean(res))/sd(res)
kf <- suppressMessages(forecast(auto.arima(tsclean(k)),h=h)$mean)
gf <- suppressMessages(forecast(auto.arima(tsclean(g),stationary=TRUE),h=h)$mean); gg <- c(g,gf)
Mf <- array(NA,c(h,nc)); Mfs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- exp(a[j]+kf[i]+gg[i-j+nc+nr]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- M[nr,j]*exp(kf[i]-k[nr]+gg[i-j+nc+nr]-g[nr-j+nc]) }}}
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M=M,h=h,jumpoff=jumpoff,alpha=a,kappa=k,gamma=g,standardresiduals=res,forecast=Mf,smoothforecast=Mfs),
class="APCS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.APCS <- function(object,...) {
list(alpha=object$alpha,kappa=object$kappa,gamma=object$gamma)
}

#' @export
forecast.APCS <- function(object,...) {
object$smoothforecast
}

#' @export
plot.APCS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(3,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha,xlab="age",ylab="alpha",pch=16,cex=0.5,bty="n")
plot(c(1:nrow(x$M)),x$kappa,xlab="year",ylab="kappa",pch=16,cex=0.5,bty="n")
plot(c(1:(nrow(x$M)+ncol(x$M)-1)),x$gamma,xlab="cohort",ylab="gamma",pch=16,cex=0.5,bty="n")
colband <- colorRampPalette(c("black","grey","white"))
image(x=x$x,y=c(1:nrow(x$M)),z=t(x$standardresiduals),col=colband(6),xlab="age",ylab="year",main="standardised residuals",cex.main=1,font.main=1)
plot(x$x,log(x$M[1,]),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M),log(x$smoothforecast)),max(log(x$M),log(x$smoothforecast))))
points(x$x,log(x$M[nrow(x$M),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
residuals.APCS <- function(object,...) {
object$standardresiduals
}

#' @export
simulate.APCS <- function(object,nsim=10,seed=123,...) {
if (!is.numeric(nsim)||nsim!=floor(nsim)||nsim<1||!is.numeric(seed)||seed!=floor(seed)||seed<1) { stop("nsim and seed must be positive integers") }
M <- object$M; h <- object$h; jumpoff <- object$jumpoff
a <- object$alpha; k <- object$kappa; g <- object$gamma
nr <- nrow(M); nc <- ncol(M)
amat <- matrix(a,nr,nc,byrow=TRUE)
kmat <- matrix(k,nr,nc,byrow=FALSE)
gmat <- g[row(M)-col(M)+nc]
res <- log(M)-amat-kmat-gmat
s <- sd(res)
f1 <- suppressMessages(auto.arima(tsclean(k)))
f2 <- suppressMessages(auto.arima(tsclean(g),stationary=TRUE))
Msim <- array(NA,c(h,nc,nsim))
set.seed(seed)
if (jumpoff==1) {
for (z in 1:nsim) { 
ksim <- simulate(f1,nsim=h)
gsim <- simulate(f2,nsim=h)
ggsim <- c(g,gsim)
for (i in 1:h) { for (j in 1:nc) { Msim[i,j,z] <- exp(a[j]+ksim[i]+ggsim[i-j+nc+nr]+rnorm(1,0,s)) }}
}} 
if (jumpoff==2) {
for (z in 1:nsim) {
ksim <- simulate(f1,nsim=h)
gsim <- simulate(f2,nsim=h)
ggsim <- c(g,gsim)
for (i in 1:h) { for (j in 1:nc) { Msim[i,j,z] <- M[nr,j]*exp(ksim[i]-k[nr]+ggsim[i-j+nc+nr]-g[nr-j+nc]+rnorm(1,0,s)) }}
}}
list(Msim=Msim,sigma=s)
}
