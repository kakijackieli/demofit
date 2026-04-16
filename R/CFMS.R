#' Common factor model
#'
#' Fits and forecasts mortality rates of two populations using common factor model.
#'
#' @param x vector of ages.
#' @param M1 matrix of mortality rates of population 1 (rows as years and columns as ages).
#' @param M2 matrix of mortality rates of population 2 (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed mortality rates in final year (default = 1). 
#'
#' @details
#' The common factor model (CFM) is specified as 
#' 
#' \eqn{ln(m_{x,t,i}) = \alpha_{x,i} + B_x K_t + \beta_{x,i} \kappa_{t,i} + \epsilon_{x,t,i}}.
#'
#' The model is estimated by Newton updating scheme and is forecasted by ARIMA applied to \eqn{K_t} and \eqn{\kappa_{t,i}}. Constraints include sum of \eqn{B_x} is one, sum of \eqn{K_t} is zero, sum of \eqn{\beta_{x,i}} is one, and sum of \eqn{\kappa_{t,i}} is zero. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom stats fitted prcomp sd
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class CFMS with associated S3 methods coef, forecast, plot (which = 1 gives parameter estimates; which = 2 gives residuals and forecasts), and residuals.
#'
#' @references
#' Li, N. and Lee, R. (2005). Coherent mortality forecasts for a group of populations: An extension of the Lee-Carter method. Demography, 42(3), 575-594.
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
#' b1 <- c(0.0012,-0.0033,0.0523,0.0161,0.0529,0.0220,0.0312,0.0437,0.0709,0.0444,
#' 0.0398,0.0361,0.0403,0.0396,0.0506,0.0315,0.0428,0.0261,0.0384,0.0388,
#' 0.0300,0.0269,0.0275,0.0256,0.0239,0.0421,0.0314,0.0284,0.0174,0.0314)
#' k1 <- c(-1.24,-1.38,-3.48,-2.51,-1.32,-1.90,-3.42,-0.94,0.24,-0.48,
#' -0.26,2.70,1.39,-0.46,1.74,2.53,0.90,1.43,0.76,2.48,
#' 0.74,2.32,0.42,1.69,-0.64,1.30,0.19,-0.69,-1.11,-1.01)
#' b2 <- c(-0.0014,0.0272,0.0083,0.0273,0.0209,0.0253,0.0144,0.0333,0.0460,0.0439,
#' 0.0439,0.0674,0.0331,0.0443,0.0312,0.0240,0.0570,0.0312,0.0403,0.0376,
#' 0.0500,0.0289,0.0466,0.0418,0.0349,0.0149,0.0366,0.0178,0.0361,0.0372)
#' k2 <- c(2.35,0.62,-0.38,0.12,0.00,0.80,-1.39,0.38,2.47,0.40,
#' 0.76,3.06,1.42,-0.73,0.79,1.94,0.12,0.60,-0.43,0.29,
#' 0.17,0.98,-1.01,-0.13,-2.46,-1.24,-1.65,-2.48,-2.32,-3.06)
#' set.seed(123)
#' M1 <- exp(outer(k1,b1)+outer(K,B)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' M2 <- exp(outer(k2,b2)+outer(K,B)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
#' fit <- CFMS(x=x,M1=M1,M2=M2,curve="makeham",h=30,jumpoff=2)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
CFMS <- function(x,M1,M2,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10,jumpoff=1) {
if (!is.numeric(x)||!is.numeric(M1)||!is.numeric(M2)) { stop("x and M1 and M2 must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M1)||!is.matrix(M2)) stop("M1 and M2 must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M1)||length(x)!=ncol(M2)) stop("the number of ages must match the number of columns of M1 and M2")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(M1<=0)||any(M2<=0)) { stop("all M1 and M2 values must be positive") }
if (nrow(M1)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
nr <- nrow(M1); nc <- ncol(M1)
PCA <- prcomp(log((M1+M2)/2),center=TRUE,scale=FALSE)
a1 <- numeric(); for (j in 1:nc) { a1[j] <- mean(log(M1[,j])) }
a2 <- numeric(); for (j in 1:nc) { a2[j] <- mean(log(M2[,j])) }
B <- PCA$rotation[,1]/sum(PCA$rotation[,1])
K <- PCA$x[,1]*sum(PCA$rotation[,1])
PCA1 <- prcomp(log(M1)-matrix(B,nr,nc,byrow=TRUE)*matrix(K,nr,nc,byrow=FALSE),center=TRUE,scale=FALSE)
b1 <- PCA1$rotation[,1]/sum(PCA1$rotation[,1])
k1 <- PCA1$x[,1]*sum(PCA1$rotation[,1])
PCA2 <- prcomp(log(M2)-matrix(B,nr,nc,byrow=TRUE)*matrix(K,nr,nc,byrow=FALSE),center=TRUE,scale=FALSE)
b2 <- PCA2$rotation[,1]/sum(PCA2$rotation[,1])
k2 <- PCA2$x[,1]*sum(PCA2$rotation[,1])
olde <- 1000000; tol <- 1e-8
for (z in 1:200) {
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
a1 <- a1+colSums(log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)/nr
a2 <- a2+colSums(log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)/nr
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
B <- B+(colSums((log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)*Kmat)+colSums((log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)*Kmat))/sum(K^2)/2; B <- B/sum(B)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
K <- K+(rowSums((log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)*Bmat)+rowSums((log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)*Bmat))/sum(B^2)/2; K <- K-mean(K)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
b1 <- b1+colSums((log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)*k1mat)/sum(k1^2); b1 <- b1/sum(b1)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
k1 <- k1+rowSums((log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)*b1mat)/sum(b1^2); k1 <- k1-mean(k1)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
b2 <- b2+colSums((log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)*k2mat)/sum(k2^2); b2 <- b2/sum(b2)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
k2 <- k2+rowSums((log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)*b2mat)/sum(b2^2); k2 <- k2-mean(k2)
a1mat <- matrix(a1,nr,nc,byrow=TRUE); a2mat <- matrix(a2,nr,nc,byrow=TRUE); Bmat <- matrix(B,nr,nc,byrow=TRUE); Kmat <- matrix(K,nr,nc,byrow=FALSE); b1mat <- matrix(b1,nr,nc,byrow=TRUE); k1mat <- matrix(k1,nr,nc,byrow=FALSE); b2mat <- matrix(b2,nr,nc,byrow=TRUE); k2mat <- matrix(k2,nr,nc,byrow=FALSE)
newe <- sum((log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat)^2)+sum((log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat)^2)
if (z>10&&olde-newe<tol&&olde>newe) { break }
olde <- newe
}
res1 <- log(M1)-a1mat-Bmat*Kmat-b1mat*k1mat
res2 <- log(M2)-a2mat-Bmat*Kmat-b2mat*k2mat
res1 <- (res1-mean(res1))/sd(res1)
res2 <- (res2-mean(res2))/sd(res2)
Kf <- suppressMessages(forecast(auto.arima(tsclean(K)),h=h)$mean)
k1f <- suppressMessages(forecast(auto.arima(tsclean(k1),stationary=TRUE),h=h)$mean)
k2f <- suppressMessages(forecast(auto.arima(tsclean(k2),stationary=TRUE),h=h)$mean)
M1f <- array(NA,c(h,nc)); M1fs <- array(NA,c(h,nc))
M2f <- array(NA,c(h,nc)); M2fs <- array(NA,c(h,nc))
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- exp(a1[j]+B[j]*Kf[i]+b1[j]*k1f[i]) }}}
if (jumpoff==1) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- exp(a2[j]+B[j]*Kf[i]+b2[j]*k2f[i]) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M1f[i,j] <- M1[nr,j]*exp(B[j]*(Kf[i]-K[nr])+b1[j]*(k1f[i]-k1[nr])) }}}
if (jumpoff==2) { for (i in 1:h) { for (j in 1:nc) { M2f[i,j] <- M2[nr,j]*exp(B[j]*(Kf[i]-K[nr])+b2[j]*(k2f[i]-k2[nr])) }}}
for (i in 1:h) { M1fs[i,] <- fitted(MC(x=x,m=M1f[i,],curve=curve)) }
for (i in 1:h) { M2fs[i,] <- fitted(MC(x=x,m=M2f[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M1=M1,M2=M2,h=h,jumpoff=jumpoff,alpha1=a1,alpha2=a2,B=B,K=K,beta1=b1,kappa1=k1,beta2=b2,kappa2=k2,standardresiduals1=res1,standardresiduals2=res2,forecast1=M1f,forecast2=M2f,smoothforecast1=M1fs,smoothforecast2=M2fs),
class="CFMS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.CFMS <- function(object,...) {
list(alpha1=object$alpha1,alpha2=object$alpha2,B=object$B,K=object$K,beta1=object$beta1,kappa1=object$kappa1,beta2=object$beta2,kappa2=object$kappa2)
}

#' @export
forecast.CFMS <- function(object,...) {
list(smoothforecast1=object$smoothforecast1,smoothforecast2=object$smoothforecast2)
}

#' @export
plot.CFMS <- function(x,which=1,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) {
par(mfrow=c(4,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
plot(x$x,x$alpha1,xlab="age",ylab="alpha1",pch=16,cex=0.5,bty="n")
plot(x$x,x$alpha2,xlab="age",ylab="alpha2",pch=16,cex=0.5,bty="n")
plot(x$x,x$B,xlab="age",ylab="B",pch=16,cex=0.5,bty="n")
plot(c(1:nrow(x$M1)),x$K,xlab="year",ylab="K",pch=16,cex=0.5,bty="n")
plot(x$x,x$beta1,xlab="age",ylab="beta1",pch=16,cex=0.5,bty="n")
plot(c(1:nrow(x$M1)),x$kappa1,xlab="year",ylab="kappa1",pch=16,cex=0.5,bty="n")
plot(x$x,x$beta2,xlab="age",ylab="beta2",pch=16,cex=0.5,bty="n")
plot(c(1:nrow(x$M2)),x$kappa2,xlab="year",ylab="kappa2",pch=16,cex=0.5,bty="n")
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
residuals.CFMS <- function(object,...) {
list(standardresiduals1=object$standardresiduals1,standardresiduals2=object$standardresiduals2)
}
