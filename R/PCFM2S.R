#' Poisson common factor model with global age pattern
#'
#' Fits and forecasts mortality rates of two populations using common factor model with global age pattern with Poisson assumption.
#'
#' @param x vector of ages.
#' @param D1 matrix of death counts of population 1 (rows as years and columns as ages).
#' @param D2 matrix of death counts of population 2 (rows as years and columns as ages).
#' @param E1 matrix of mid-year exposures of population 1 (rows as years and columns as ages).
#' @param E2 matrix of mid-year exposures of population 2 (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The common factor model with global age pattern with Poisson assumption is specified as 
#' 
#' \eqn{ln(m_{x,t,i}) = \alpha_{x,i} + B_x K_t + \beta_x \kappa_{t,i}} and \eqn{D_{x,t,i} ~ Poisson(E_{x,t,i} m_{x,t,i})}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{K_t} and \eqn{\kappa_{t,i}}. Constraints include sum of \eqn{B_x} is one, sum of \eqn{K_t} is zero, sum of \eqn{\beta_x} is one, and sum of \eqn{\kappa_{t,i}} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class PCFM2S with associated S3 methods coef, forecast (which = 1 for smoothed (default); which = 2 for raw), plot (which = 1 gives parameter estimates (default); which = 2 gives residuals and forecasts), and residuals.
#'
#' @references
#' Li, J., Wang, M., Liu, J., and Leung, J.W.Y. (2026). Financial valuation of retirement village via stochastic modelling of disability prevalence rates. ASTIN Bulletin, 56(2), 447-473.
#'
#' @examples
#' x <- 60:89
#' a1 <- c(-5.18,-5.12,-4.98,-4.92,-4.82,-4.73,-4.66,-4.53,-4.45,-4.35,
#' -4.26,-4.17,-4.05,-3.95,-3.84,-3.73,-3.65,-3.52,-3.40,-3.29,
#' -3.14,-3.02,-2.88,-2.76,-2.64,-2.49,-2.37,-2.25,-2.12,-2.00)
#' a2 <- c(-4.78,-4.68,-4.57,-4.49,-4.39,-4.29,-4.19,-4.10,-4.00,-3.89,
#' -3.80,-3.69,-3.60,-3.49,-3.39,-3.29,-3.17,-3.07,-2.96,-2.85,
#' -2.71,-2.62,-2.49,-2.37,-2.26,-2.14,-2.04,-1.91,-1.82,-1.72)
#' B <- c(0.0381,0.0340,0.0420,0.0389,0.0423,0.0414,0.0406,0.0393,0.0415,0.0400,
#' 0.0411,0.0362,0.0387,0.0381,0.0384,0.0385,0.0356,0.0314,0.0317,0.0337,
#' 0.0316,0.0298,0.0284,0.0270,0.0248,0.0262,0.0205,0.0215,0.0142,0.0145)
#' K <- c(9.66,9.89,10.66,9.83,9.52,7.39,7.64,6.36,2.32,4.18,
#' 2.91,-0.61,0.28,-0.38,-1.79,-3.34,-1.74,-3.50,-4.28,-4.77,
#' -4.98,-7.13,-5.09,-6.41,-5.56,-5.65,-6.12,-5.64,-7.35,-6.28)
#' b <- c(-0.00010,0.01195,0.03030,0.02170,0.03690,0.02365,0.02280,0.03850,0.05845,0.04415,
#' 0.04185,0.05175,0.03670,0.04195,0.04090,0.02775,0.04990,0.02865,0.03935,0.03820,
#' 0.04000,0.02790,0.03705,0.03370,0.02940,0.02850,0.03400,0.02310,0.02675,0.03430)
#' k1 <- c(-1.24,-1.38,-3.48,-2.51,-1.32,-1.90,-3.42,-0.94,0.24,-0.48,
#' -0.26,2.70,1.39,-0.46,1.74,2.53,0.90,1.43,0.76,2.48,
#' 0.74,2.32,0.42,1.69,-0.64,1.30,0.19,-0.69,-1.11,-1.01)
#' k2 <- c(2.35,0.62,-0.38,0.12,0.00,0.80,-1.39,0.38,2.47,0.40,
#' 0.76,3.06,1.42,-0.73,0.79,1.94,0.12,0.60,-0.43,0.29,
#' 0.17,0.98,-1.01,-0.13,-2.46,-1.24,-1.65,-2.48,-2.32,-3.06)
#' set.seed(123)
#' M1 <- exp(outer(k1,b)+outer(K,B)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' M2 <- exp(outer(k2,b)+outer(K,B)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' E1 <- matrix(c(107788,108036,107481,106552,104608,100104,95803,91345,84980,79557,
#' 75146,70559,65972,60898,55623,50522,47430,45895,41443,34774,
#' 30531,27754,25105,22271,19437,16888,14458,12146,10038,7994),30,30,byrow=TRUE)
#' E2 <- E1
#' D1 <- round(E1*M1)
#' D2 <- round(E2*M2)
#' fit <- PCFM2S(x=x,D1=D1,D2=D2,E1=E1,E2=E2,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
PCFM2S <- function(x,D1,D2,E1,E2,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(D1)||!is.numeric(D2)||!is.numeric(E1)||!is.numeric(E2)) { stop("x, D1, D2, E1, and E2 must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(D1)||!is.matrix(D2)||!is.matrix(E1)||!is.matrix(E2)) stop("D1, D2, E1, and E2 must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(D1)||length(x)!=ncol(D2)||length(x)!=ncol(E1)||length(x)!=ncol(E2)) stop("the number of ages must match the number of columns of D1, D2, E1, and E2")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(D1<=0)||any(D2<=0)||any(E1<=0)||any(E2<=0)) { stop("all D1, D2, E1, and E2 values must be positive") }
if (nrow(D1)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
M1 <- D1/E1; M2 <- D2/E2
nr <- nrow(M1); nc <- ncol(M1)
fit <- CFM2S(x,M1,M2,curve,h,jumpoff)
a1 <- coef(fit)$alpha1
a2 <- coef(fit)$alpha2
B <- coef(fit)$B
K <- coef(fit)$K
b <- coef(fit)$beta
k1 <- coef(fit)$kappa1
k2 <- coef(fit)$kappa2
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
a1 <- a1+colSums(D1-E1mat)/colSums(E1mat)
a2 <- a2+colSums(D2-E2mat)/colSums(E2mat)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
B <- B+(colSums((D1-E1mat)*Kmat)+colSums((D2-E2mat)*Kmat))/(colSums(E1mat*Kmat^2)+colSums(E2mat*Kmat^2)); B <- B/sum(B)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
K <- K+(rowSums((D1-E1mat)*Bmat)+rowSums((D2-E2mat)*Bmat))/(rowSums(E1mat*Bmat^2)+rowSums(E2mat*Bmat^2)); K <- K-mean(K)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
b <- b+(colSums((D1-E1mat)*k1mat)+colSums((D2-E2mat)*k2mat))/(colSums(E1mat*k1mat^2)+colSums(E2mat*k2mat^2)); b <- b/sum(b)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
k1 <- k1+rowSums((D1-E1mat)*bmat)/rowSums(E1mat*bmat^2); k1 <- k1-mean(k1)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
k2 <- k2+rowSums((D2-E2mat)*bmat)/rowSums(E2mat*bmat^2); k2 <- k2-mean(k2)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+Bmat*Kmat+bmat*k1mat); E2mat <- E2*exp(a2mat+Bmat*Kmat+bmat*k2mat)
newe <- sum(D1*log(D1/E1mat)-D1+E1mat)+sum(D2*log(D2/E2mat)-D2+E2mat)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
dis1 <- sum(2*(D1*log(D1/E1mat)-D1+E1mat))/(nr*nc-nc-nc-nr*1.5+2.5)
res1 <- sign(D1-E1mat)*sqrt(2*(D1*log(D1/E1mat)-D1+E1mat)/dis1)
dis2 <- sum(2*(D2*log(D2/E2mat)-D2+E2mat))/(nr*nc-nc-nc-nr*1.5+2.5)
res2 <- sign(D2-E2mat)*sqrt(2*(D2*log(D2/E2mat)-D2+E2mat)/dis2)
Kf <- suppressMessages(forecast(auto.arima(tsclean(K)),h=h)$mean)
k1f <- suppressMessages(forecast(auto.arima(tsclean(k1),stationary=TRUE),h=h)$mean)
k2f <- suppressMessages(forecast(auto.arima(tsclean(k2),stationary=TRUE),h=h)$mean)
M1f <- array(NA,c(h,nc)); M1fs <- array(NA,c(h,nc))
M2f <- array(NA,c(h,nc)); M2fs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- exp(a1[j]+B[j]*Kf[i]+b[j]*k1f[i]) }}}
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- exp(a2[j]+B[j]*Kf[i]+b[j]*k2f[i]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- M1[nr,j]*exp(B[j]*(Kf[i]-K[nr])+b[j]*(k1f[i]-k1[nr])) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- M2[nr,j]*exp(B[j]*(Kf[i]-K[nr])+b[j]*(k2f[i]-k2[nr])) }}}
for (i in 1:h) { M1fs[i,] <- fitted(MC(x=x,m=M1f[i,],curve=curve)) }
for (i in 1:h) { M2fs[i,] <- fitted(MC(x=x,m=M2f[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,D1=D1,D2=D2,E1=E1,E2=E2,M1=M1,M2=M2,h=h,jumpoff=jumpoff,alpha1=a1,alpha2=a2,B=B,K=K,beta=b,kappa1=k1,kappa2=k2,standardresiduals1=res1,standardresiduals2=res2,dispersion1=dis1,dispersion2=dis2,forecast1=M1f,forecast2=M2f,smoothforecast1=M1fs,smoothforecast2=M2fs),
class="PCFM2S"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.PCFM2S <- function(object,...) {
list(alpha1=object$alpha1,alpha2=object$alpha2,B=object$B,K=object$K,beta=object$beta,kappa1=object$kappa1,kappa2=object$kappa2)
}

#' @export
forecast.PCFM2S <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { list(smoothforecast1=object$smoothforecast1,smoothforecast2=object$smoothforecast2) }
else if (which==2) { list(forecast1=object$forecast1,forecast2=object$forecast2) }
}

#' @export
plot.PCFM2S <- function(x,which=1,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) {
par(mfrow=c(4,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha1,xlab="age",ylab="alpha1",bty="n",type="l")
plot(x$x,x$alpha2,xlab="age",ylab="alpha2",bty="n",type="l")
plot(x$x,x$B,xlab="age",ylab="B",bty="n",type="l")
plot(c(1:nrow(x$M1)),x$K,xlab="year",ylab="K",bty="n",type="l")
plot(x$x,x$beta,xlab="age",ylab="beta",bty="n",type="l")
plot(c(1:nrow(x$M1)),x$kappa1,xlab="year",ylab="kappa1",bty="n",type="l")
plot(c(1:nrow(x$M2)),x$kappa2,xlab="year",ylab="kappa2",bty="n",type="l")
}
if (which==2) {
par(mfrow=c(2,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
colband <- colorRampPalette(c("black","grey","white"))
image(x=x$x,y=c(1:nrow(x$M1)),z=t(x$standardresiduals1),col=colband(6),xlab="age",ylab="year",main="standardised residuals 1",cex.main=1,font.main=1)
image(x=x$x,y=c(1:nrow(x$M2)),z=t(x$standardresiduals2),col=colband(6),xlab="age",ylab="year",main="standardised residuals 2",cex.main=1,font.main=1)
plot(x$x,log(x$M1[1,]),xlab="age",ylab="log death rate 1",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M1),log(x$smoothforecast1)),max(log(x$M1),log(x$smoothforecast1))))
points(x$x,log(x$M1[nrow(x$M1),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast1[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
plot(x$x,log(x$M2[1,]),xlab="age",ylab="log death rate 2",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M2),log(x$smoothforecast2)),max(log(x$M2),log(x$smoothforecast2))))
points(x$x,log(x$M2[nrow(x$M2),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast2[x$h,]))
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}}

#' @export
residuals.PCFM2S <- function(object,...) {
list(standardresiduals1=object$standardresiduals1,standardresiduals2=object$standardresiduals2)
}
