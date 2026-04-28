#' CBD with cohort model
#'
#' Fits and forecasts mortality rates using CBD with cohort model.
#'
#' @param x vector of ages.
#' @param M matrix of mortality rates (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, perks, weibull, beard, martinelle, thatcher, gompertz2, makeham2, perks2, weibull2, beard2, martinelle2, thatcher2, where first 7 curves' parameters are unconstrained and last 7 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The CBD with cohort (M6) model is specified as 
#' 
#' \eqn{ln(m_{x,t}) = \kappa_{1,t} + \kappa_{2,t} (x-\bar{x}) + \gamma_{t-x} + \epsilon_{x,t}}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_{1,t}}, \eqn{\kappa_{2,t}}, and \eqn{\gamma_c}. Constraints include sum of \eqn{\gamma_c} is zero and sum of \eqn{c\gamma_c} is zero. It is designed for ages 50-90.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted lm sd simulate rnorm
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class CBDCS with associated S3 methods coef, forecast, plot, residuals, and simulate (nsim for setting number of simulations; seed for initialising random number generator).
#'
#' @references
#' Cairns, A.J.G., Blake, D., Dowd, K., Coughlan, G.D., Epstein, D., Ong, A., and Balevich, I. (2009). A quantitative comparison of stochastic mortality models using data from England and Wales and the United States. North American Actuarial Journal, 13(1), 1-35.
#'
#' @examples
#' x <- 60:89
#' k1 <- -2.97-0.0245*(0:29)
#' k2 <- 0.101+0.000345*(0:29)
#' set.seed(123)
#' M <- exp(matrix(k1,nrow=30,ncol=30,byrow=FALSE)+outer(k2,(x-mean(x)))+rnorm(900,0,0.035))
#' fit <- CBDCS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
CBDCS <- function(x,M,curve=c("gompertz","makeham","perks","weibull","beard","martinelle","thatcher","gompertz2","makeham2","perks2","weibull2","beard2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(M)) { stop("x and M must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M)) stop("M must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M)) stop("the number of ages must match the number of columns of M")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if ((min(x)<50)||(max(x)>90)) { stop("x must be in between 50 and 90") }
if (any(M<=0)) { stop("all M values must be positive") }
if (nrow(M)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
nr <- nrow(M); nc <- ncol(M)
k1 <- numeric(); k2 <- numeric()
xx <- x-mean(x); for (i in 1:nr) { yy <- log(M[i,]); k1[i] <- as.numeric(lm(yy~xx)$coef[1]); k2[i] <- as.numeric(lm(yy~xx)$coef[2]) }
g <- rep(0,nr+nc-1)
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]; xm <- matrix(x-mean(x),nr,nc,byrow=TRUE)
k1 <- k1+rowSums(log(M)-k1mat-k2mat*xm-gmat)/nc
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]; xm <- matrix(x-mean(x),nr,nc,byrow=TRUE)
k2 <- k2+rowSums((log(M)-k1mat-k2mat*xm-gmat)*xm)/sum((x-mean(x))^2)
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]; xm <- matrix(x-mean(x),nr,nc,byrow=TRUE)
g <- g+as.numeric(tapply(as.vector(log(M)-k1mat-k2mat*xm-gmat),as.vector(row(M)-col(M)+nc),mean)); g <- g-mean(g)
cc <- c(1:(nr+nc-1)); C <- cbind(1,cc); g <- as.numeric(g-C%*%qr.solve(t(C)%*%C)%*%t(C)%*%g)
k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); gmat <- g[row(M)-col(M)+nc]; xm <- matrix(x-mean(x),nr,nc,byrow=TRUE)
newe <- sum((log(M)-k1mat-k2mat*xm-gmat)^2)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
res <- array(NA,c(nr,nc))
for (i in 1:nr) { for (j in 1:nc) { res[i,j] <- log(M[i,j])-k1[i]-k2[i]*(x[j]-mean(x))-g[i-j+nc] }}
res <- (res-mean(res))/sd(res)
k1f <- suppressMessages(forecast(auto.arima(tsclean(k1)),h=h)$mean)
k2f <- suppressMessages(forecast(auto.arima(tsclean(k2)),h=h)$mean)
gf <- suppressMessages(forecast(auto.arima(tsclean(g),stationary=TRUE),h=h)$mean); gg <- c(g,gf)
Mf <- array(NA,c(h,nc)); Mfs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- exp(k1f[i]+k2f[i]*(x[j]-mean(x))+gg[i-j+nc+nr]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { Mf[i,j] <- M[nr,j]*exp(k1f[i]-k1[nr]+(k2f[i]-k2[nr])*(x[j]-mean(x))+gg[i-j+nc+nr]-g[nr-j+nc]) }}}
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M=M,h=h,jumpoff=jumpoff,kappa1=k1,kappa2=k2,gamma=g,standardresiduals=res,forecast=Mf,smoothforecast=Mfs),
class="CBDCS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.CBDCS <- function(object,...) {
list(kappa1=object$kappa1,kappa2=object$kappa2,gamma=object$gamma)
}

#' @export
forecast.CBDCS <- function(object,...) {
object$smoothforecast
}

#' @export
plot.CBDCS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(3,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(c(1:nrow(x$M)),x$kappa1,xlab="year",ylab="kappa1",pch=16,cex=0.5,bty="n")
plot(c(1:nrow(x$M)),x$kappa2,xlab="year",ylab="kappa2",pch=16,cex=0.5,bty="n")
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
residuals.CBDCS <- function(object,...) {
object$standardresiduals
}

#' @export
simulate.CBDCS <- function(object,nsim=10,seed=123,...) {
if (!is.numeric(nsim)||nsim!=floor(nsim)||nsim<1||!is.numeric(seed)||seed!=floor(seed)||seed<1) { stop("nsim and seed must be positive integers") }
x <- object$x; M <- object$M; h <- object$h; jumpoff <- object$jumpoff
k1 <- object$kappa1; k2 <- object$kappa2; g <- object$gamma
nr <- nrow(M); nc <- ncol(M)
k1mat <- matrix(k1,nr,nc,byrow=FALSE)
k2mat <- matrix(k2,nr,nc,byrow=FALSE)
gmat <- g[row(M)-col(M)+nc]
xmat <- matrix(x-mean(x),nr,nc,byrow=TRUE)
res <- log(M)-k1mat-k2mat*xmat-gmat
s <- sd(res)
f1 <- suppressMessages(auto.arima(tsclean(k1)))
f2 <- suppressMessages(auto.arima(tsclean(k2)))
f3 <- suppressMessages(auto.arima(tsclean(g),stationary=TRUE))
Msim <- array(NA,c(h,nc,nsim))
set.seed(seed)
if (jumpoff==1) {
for (z in 1:nsim) { 
k1sim <- simulate(f1,nsim=h)
k2sim <- simulate(f2,nsim=h)
gsim <- simulate(f3,nsim=h)
ggsim <- c(g,gsim)
for (i in 1:h) { for (j in 1:nc) { Msim[i,j,z] <- exp(k1sim[i]+k2sim[i]*(x[j]-mean(x))+ggsim[i-j+nc+nr]+rnorm(1,0,s)) }}
}} 
if (jumpoff==2) {
for (z in 1:nsim) {
k1sim <- simulate(f1,nsim=h)
k2sim <- simulate(f2,nsim=h)
gsim <- simulate(f3,nsim=h)
ggsim <- c(g,gsim)
for (i in 1:h) { for (j in 1:nc) { Msim[i,j,z] <- M[nr,j]*exp(k1sim[i]-k1[nr]+(k2sim[i]-k2[nr])*(x[j]-mean(x))+ggsim[i-j+nc+nr]-g[nr-j+nc]+rnorm(1,0,s)) }}
}}
list(Msim=Msim,sigma=s)
}
