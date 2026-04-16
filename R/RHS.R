#' Renshaw-Haberman model
#'
#' Fits and forecasts mortality rates using Renshaw-Haberman model.
#'
#' @param x vector of ages.
#' @param M matrix of mortality rates (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The Renshaw-Haberman (RH) model is specified as 
#' 
#' \eqn{ln(m_{x,t}) = \alpha_x + \beta_x \kappa_t + \gamma_{t-x} + \epsilon_{x,t}}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_t} and \eqn{\gamma_c}. Constraints include sum of \eqn{\beta_x} is one, sum of \eqn{\kappa_t} is zero, and sum of \eqn{\gamma_c} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class RHS with associated S3 methods coef, forecast, plot, and residuals.
#'
#' @references
#' Haberman, S. and Renshaw, A. (2011). A comparative study of parametric mortality projection models. Insurance: Mathematics and Economics, 48(1), 35-55.
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
#' fit <- RHS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
RHS <- function(x,M,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
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
PCA <- prcomp(log(M),center=TRUE,scale=FALSE)
a <- numeric(); for (j in 1:ncol(M)) { a[j] <- mean(log(M[,j])) }
b <- PCA$rotation[,1]/sum(PCA$rotation[,1])
k <- PCA$x[,1]*sum(PCA$rotation[,1])
g <- rep(0,nrow(M)+ncol(M)-1)
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
amat <- matrix(a,nrow(M),ncol(M),byrow=TRUE); bmat <- matrix(b,nrow(M),ncol(M),byrow=TRUE); kmat <- matrix(k,nrow(M),ncol(M),byrow=FALSE); gmat <- g[row(M)-col(M)+ncol(M)]
a <- a+colSums(log(M)-amat-bmat*kmat-gmat)/nrow(M)
amat <- matrix(a,nrow(M),ncol(M),byrow=TRUE); bmat <- matrix(b,nrow(M),ncol(M),byrow=TRUE); kmat <- matrix(k,nrow(M),ncol(M),byrow=FALSE); gmat <- g[row(M)-col(M)+ncol(M)]
b <- b+colSums((log(M)-amat-bmat*kmat-gmat)*kmat)/sum(k^2); b <- b/sum(b)
amat <- matrix(a,nrow(M),ncol(M),byrow=TRUE); bmat <- matrix(b,nrow(M),ncol(M),byrow=TRUE); kmat <- matrix(k,nrow(M),ncol(M),byrow=FALSE); gmat <- g[row(M)-col(M)+ncol(M)]
k <- k+rowSums((log(M)-amat-bmat*kmat-gmat)*bmat)/sum(b^2); k <- k-mean(k)
amat <- matrix(a,nrow(M),ncol(M),byrow=TRUE); bmat <- matrix(b,nrow(M),ncol(M),byrow=TRUE); kmat <- matrix(k,nrow(M),ncol(M),byrow=FALSE); gmat <- g[row(M)-col(M)+ncol(M)]
g <- g+as.numeric(tapply(as.vector(log(M)-amat-bmat*kmat-gmat),as.vector(row(M)-col(M)+ncol(M)),mean)); g <- g-mean(g)
amat <- matrix(a,nrow(M),ncol(M),byrow=TRUE); bmat <- matrix(b,nrow(M),ncol(M),byrow=TRUE); kmat <- matrix(k,nrow(M),ncol(M),byrow=FALSE); gmat <- g[row(M)-col(M)+ncol(M)]
newe <- sum((log(M)-amat-bmat*kmat-gmat)^2)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
res <- array(NA,c(nrow(M),ncol(M)))
for (i in 1:nrow(M)) { for (j in 1:ncol(M)) { res[i,j] <- log(M[i,j])-a[j]-b[j]*k[i]-g[i-j+ncol(M)] }}
res <- (res-mean(res))/sd(res)
kf <- suppressMessages(forecast(auto.arima(tsclean(k)),h=h)$mean)
gf <- suppressMessages(forecast(auto.arima(tsclean(g),stationary=TRUE),h=h)$mean); gg <- c(g,gf)
Mf <- array(NA,c(h,ncol(M))); Mfs <- array(NA,c(h,ncol(M)))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:ncol(M)) { Mf[i,j] <- exp(a[j]+b[j]*kf[i]+gg[i-j+ncol(M)+nrow(M)]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:ncol(M)) { Mf[i,j] <- M[nrow(M),j]*exp(b[j]*(kf[i]-k[nrow(M)])+gg[i-j+ncol(M)+nrow(M)]-g[nrow(M)-j+ncol(M)]) }}}
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M=M,h=h,jumpoff=jumpoff,alpha=a,beta=b,kappa=k,gamma=g,standardresiduals=res,forecast=Mf,smoothforecast=Mfs),
class="RHS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.RHS <- function(object,...) {
list(alpha=object$alpha,beta=object$beta,kappa=object$kappa,gamma=object$gamma)
}

#' @export
forecast.RHS <- function(object,...) {
object$smoothforecast
}

#' @export
plot.RHS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(3,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha,xlab="age",ylab="alpha",pch=16,cex=0.5,bty="n")
plot(x$x,x$beta,xlab="age",ylab="beta",pch=16,cex=0.5,bty="n")
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
residuals.RHS <- function(object,...) {
object$standardresiduals
}
