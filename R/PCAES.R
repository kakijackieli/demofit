#' Poisson common age effect model
#'
#' Fits and forecasts mortality rates of two populations using common age effect model with Poisson assumption.
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
#' The common age effect (CAE) model with Poisson assumption is specified as 
#' 
#' \eqn{ln(m_{x,t,i}) = \alpha_{x,i} + \beta_x \kappa_{t,i}} and \eqn{D_{x,t,i} ~ Poisson(E_{x,t,i} m_{x,t,i})}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{\kappa_{t,i}}. Constraints include sum of \eqn{\beta_x} is one and sum of \eqn{\kappa_{t,i}} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class PCAES with associated S3 methods coef, forecast (which = 1 for smoothed (default); which = 2 for raw), plot (which = 1 gives parameter estimates (default); which = 2 gives residuals and forecasts), and residuals.
#'
#' @references
#' Kleinow, T. (2015). A common age effect model for the mortality of multiple populations. Insurance: Mathematics and Economics, 63(C), 147-152.
#'
#' @examples
#' x <- 60:89
#' a1 <- c(-5.18,-5.12,-4.98,-4.92,-4.82,-4.73,-4.66,-4.53,-4.45,-4.35,
#' -4.26,-4.17,-4.05,-3.95,-3.84,-3.73,-3.65,-3.52,-3.40,-3.29,
#' -3.14,-3.02,-2.88,-2.76,-2.64,-2.49,-2.37,-2.25,-2.12,-2.00)
#' a2 <- c(-4.78,-4.68,-4.57,-4.49,-4.39,-4.29,-4.19,-4.10,-4.00,-3.89,
#' -3.80,-3.69,-3.60,-3.49,-3.39,-3.29,-3.17,-3.07,-2.96,-2.85,
#' -2.71,-2.62,-2.49,-2.37,-2.26,-2.14,-2.04,-1.91,-1.82,-1.72)
#' b <- c(0.0381,0.0340,0.0420,0.0389,0.0423,0.0414,0.0406,0.0393,0.0415,0.0400,
#' 0.0411,0.0362,0.0387,0.0381,0.0384,0.0385,0.0356,0.0314,0.0317,0.0337,
#' 0.0316,0.0298,0.0284,0.0270,0.0248,0.0262,0.0205,0.0215,0.0142,0.0145)
#' k1 <- c(8.68,8.34,7.99,6.87,8.18,5.73,4.83,5.20,2.74,3.22,
#' 2.99,1.59,1.67,-0.65,-0.39,-1.07,-0.95,-2.78,-3.46,-2.45,
#' -4.12,-4.66,-4.98,-4.58,-6.30,-4.39,-5.56,-6.52,-8.26,-6.92)
#' k2 <- c(11.81,11.01,10.59,10.40,9.75,8.15,6.07,6.45,4.60,4.57,
#' 4.15,1.49,1.77,-1.08,-1.44,-0.96,-1.66,-2.25,-4.67,-4.62,
#' -4.38,-6.37,-6.27,-6.91,-8.22,-7.35,-8.39,-7.87,-9.72,-8.65)
#' set.seed(123)
#' M1 <- exp(outer(k1,b)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' M2 <- exp(outer(k2,b)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' E1 <- matrix(c(107788,108036,107481,106552,104608,100104,95803,91345,84980,79557,
#' 75146,70559,65972,60898,55623,50522,47430,45895,41443,34774,
#' 30531,27754,25105,22271,19437,16888,14458,12146,10038,7994),30,30,byrow=TRUE)
#' E2 <- E1
#' D1 <- round(E1*M1)
#' D2 <- round(E2*M2)
#' fit <- PCAES(x=x,D1=D1,D2=D2,E1=E1,E2=E2,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
PCAES <- function(x,D1,D2,E1,E2,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
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
fit <- CAES(x,M1,M2,curve,h,jumpoff)
a1 <- coef(fit)$alpha1
a2 <- coef(fit)$alpha2
b <- coef(fit)$beta
k1 <- coef(fit)$kappa1
k2 <- coef(fit)$kappa2
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+bmat*k1mat); E2mat <- E2*exp(a2mat+bmat*k2mat)
a1 <- a1+colSums(D1-E1mat)/colSums(E1mat)
a2 <- a2+colSums(D2-E2mat)/colSums(E2mat)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+bmat*k1mat); E2mat <- E2*exp(a2mat+bmat*k2mat)
b <- b+(colSums((D1-E1mat)*k1mat)+colSums((D2-E2mat)*k2mat))/(colSums(E1mat*k1mat^2)+colSums(E2mat*k2mat^2)); b <- b/sum(b)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+bmat*k1mat); E2mat <- E2*exp(a2mat+bmat*k2mat)
k1 <- k1+rowSums((D1-E1mat)*bmat)/rowSums(E1mat*bmat^2); k1 <- k1-mean(k1)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+bmat*k1mat); E2mat <- E2*exp(a2mat+bmat*k2mat)
k2 <- k2+rowSums((D2-E2mat)*bmat)/rowSums(E2mat*bmat^2); k2 <- k2-mean(k2)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); bmat <- matrix(b,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); k2mat <- matrix(k2,nr,nc,byrow=FALSE); E1mat <- E1*exp(a1mat+bmat*k1mat); E2mat <- E2*exp(a2mat+bmat*k2mat)
newe <- sum(D1*log(D1/E1mat)-D1+E1mat)+sum(D2*log(D2/E2mat)-D2+E2mat)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
dis1 <- sum(2*(D1*log(D1/E1mat)-D1+E1mat))/(nr*nc-nc-nc*0.5-nr+1.5)
res1 <- sign(D1-E1mat)*sqrt(2*(D1*log(D1/E1mat)-D1+E1mat)/dis1)
dis2 <- sum(2*(D2*log(D2/E2mat)-D2+E2mat))/(nr*nc-nc-nc*0.5-nr+1.5)
res2 <- sign(D2-E2mat)*sqrt(2*(D2*log(D2/E2mat)-D2+E2mat)/dis2)
k1f <- suppressMessages(forecast(auto.arima(tsclean(k1)),h=h)$mean)
k2f <- suppressMessages(forecast(auto.arima(tsclean(k2)),h=h)$mean)
M1f <- array(NA,c(h,nc)); M1fs <- array(NA,c(h,nc))
M2f <- array(NA,c(h,nc)); M2fs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- exp(a1[j]+b[j]*k1f[i]) }}}
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- exp(a2[j]+b[j]*k2f[i]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- M1[nr,j]*exp(b[j]*(k1f[i]-k1[nr])) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- M2[nr,j]*exp(b[j]*(k2f[i]-k2[nr])) }}}
for (i in 1:h) { M1fs[i,] <- fitted(MC(x=x,m=M1f[i,],curve=curve)) }
for (i in 1:h) { M2fs[i,] <- fitted(MC(x=x,m=M2f[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,D1=D1,D2=D2,E1=E1,E2=E2,M1=M1,M2=M2,h=h,jumpoff=jumpoff,alpha1=a1,alpha2=a2,beta=b,kappa1=k1,kappa2=k2,standardresiduals1=res1,standardresiduals2=res2,dispersion1=dis1,dispersion2=dis2,forecast1=M1f,forecast2=M2f,smoothforecast1=M1fs,smoothforecast2=M2fs),
class="PCAES"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.PCAES <- function(object,...) {
list(alpha1=object$alpha1,alpha2=object$alpha2,beta=object$beta,kappa1=object$kappa1,kappa2=object$kappa2)
}

#' @export
forecast.PCAES <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { list(smoothforecast1=object$smoothforecast1,smoothforecast2=object$smoothforecast2) }
else if (which==2) { list(forecast1=object$forecast1,forecast2=object$forecast2) }
}

#' @export
plot.PCAES <- function(x,which=1,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) {
par(mfrow=c(3,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha1,xlab="age",ylab="alpha1",bty="n",type="l")
plot(x$x,x$alpha2,xlab="age",ylab="alpha2",bty="n",type="l")
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
residuals.PCAES <- function(object,...) {
list(standardresiduals1=object$standardresiduals1,standardresiduals2=object$standardresiduals2)
}
