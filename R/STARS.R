#' Spatial-temporal autoregressive model
#'
#' Fits and forecasts mortality rates using spatial-temporal autoregressive model.
#'
#' @param x vector of ages.
#' @param M matrix of mortality rates (rows as years and columns as ages).
#' @param curve name of mortality curve for smoothing forecasted mortality rates (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher, gompertz2, makeham2, oppermann2, thiele2, wittsteinbumsted2, perks2, weibull2, vandermaen2, beard2, heligmanpollard2, rogersplanck2, siler2, martinelle2, thatcher2, where first 14 curves' parameters are unconstrained and last 14 curves' parameters are generally restricted to be positive).
#' @param h forecast horizon (default = 10).
#'
#' @details
#' The spatial-temporal autoregressive (STAR) model is specified as 
#' 
#' \eqn{ln(M_t) = \mu + R ln(M_{t-1}) + E_t}, 
#'
#' where \eqn{\mu} contains intercept parameters and \eqn{R} is banded matrix with non-zero entries in main diagonal and two subdiagonals below that.
#'
#' The model is estimated by constrained regression and is forecasted via its autoregressive nature. Constraints in \eqn{R} include: all non-zero entries are between zero and one and sum of each row is one. It can be applied to whole age range.
#'
#' @importFrom forecast auto.arima tsclean forecast
#' @importFrom NlcOptim solnl
#' @importFrom stats nlminb lm sd
#' @importFrom graphics par lines legend points image 
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' An object of class STARS with associated S3 methods coef, forecast, plot, and residuals.
#'
#' @references
#' Li, H. and Lu, Y. (2017). Coherent forecasting of mortality rates: A sparse vector-autoregression approach. ASTIN Bulletin, 47(2), 563-600.
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
#' fit <- STARS(x=x,M=M,curve="makeham",h=30)
#' coef(fit)
#' forecast::forecast(fit)
#' plot(fit)
#' residuals(fit)
#'
#' @export
STARS <- function(x,M,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2"),h=10) {
if (!is.numeric(x)||!is.numeric(M)) { stop("x and M must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M)) stop("M must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M)) stop("the number of ages must match the number of columns of M")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(M<=0)) { stop("all M values must be positive") }
if (nrow(M)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)) { stop("h must be numeric") }
if (h<1) { stop("h must be at least 1") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
u <- numeric(); a <- numeric(); b <- numeric(); R <- array(0,c(ncol(M),ncol(M)))
u[1] = (log(M[nrow(M),1])-log(M[1,1]))/(nrow(M)-1)
yy <- diff(log(M[,2])); xx <- log(M[-nrow(M),1])-log(M[-nrow(M),2])
fit <- lm(yy~xx)
u[2] <- as.numeric(fit$coef[1]); a[2] <- as.numeric(fit$coef[2])
if (a[2]<0||a[2]>1) {
error <- function(p) { sum((log(M[-1,2])-p[1]-p[2]*log(M[-nrow(M),1])-(1-p[2])*log(M[-nrow(M),2]))^2) }
iv <- c(u[2],0); ilower <- c(-Inf,0); iupper<- c(Inf,1)
fit <- nlminb(iv,error,lower=ilower,upper=iupper)
u[2] <- as.numeric(fit$par[1]); a[2] <- as.numeric(fit$par[2]) }
for (j in 3:ncol(M)) {
yy <- diff(log(M[,j])); xx1 <- log(M[-nrow(M),j-1])-log(M[-nrow(M),j]); xx2 <- log(M[-nrow(M),j-2])-log(M[-nrow(M),j])
fit <- lm(yy~xx1+xx2)
u[j] <- as.numeric(fit$coef[1]); a[j] <- as.numeric(fit$coef[2]); b[j] <- as.numeric(fit$coef[3])
if (a[j]<0||b[j]<0||a[j]+b[j]>1) {
error <- function(p) { sum((log(M[-1,j])-p[1]-p[3]*log(M[-nrow(M),j-2])-p[2]*log(M[-nrow(M),j-1])-(1-p[2]-p[3])*log(M[-nrow(M),j]))^2) }
iv <- c(u[j],0,0); ilower <- c(-Inf,0,0); iupper<- c(Inf,1,1)
fit <- nlminb(iv,error,lower=ilower,upper=iupper)
u[j] <- as.numeric(fit$par[1]); a[j] <- as.numeric(fit$par[2]); b[j] <- as.numeric(fit$par[3]) }
if (a[j]<0||b[j]<0||a[j]+b[j]>1) {
cons <- function(p) { list(ceq=NULL,c=p[2]+p[3]-1) }
fit <- solnl(iv,error,confun=cons,lb=ilower,ub=iupper)
u[j] <- as.numeric(fit$par[1]); a[j] <- as.numeric(fit$par[2]); b[j] <- as.numeric(fit$par[3]) }
}
R[1,1] <- 1
R[2,1]<- a[2]; R[2,2] <- 1-a[2]
for (j in 3:ncol(M)) { R[j,j-2] <- b[j]; R[j,j-1] <- a[j]; R[j,j] <-1-a[j]-b[j] }
res <- array(NA,c(nrow(M)-1,ncol(M)))
for (i in 1:(nrow(M)-1)) { res[i,] <- log(M[i+1,])-u-R%*%log(M[i,]) }
res <- (res-mean(res))/sd(res)
Mf <- array(NA,c(h,ncol(M))); Mfs <- array(NA,c(h,ncol(M)))
Mf[1,] <-  exp(u+R%*%log(M[nrow(M),]))
for (i in 2:h) { Mf[i,] <-  exp(u+R%*%log(Mf[i-1,])) }
for (i in 1:h) { Mfs[i,] <- fitted(MC(x=x,m=Mf[i,],curve=curve)) }
invisible(structure(
list(curve=curve,x=x,M=M,h=h,mu=u,R=R,standardresiduals=res,forecast=Mf,smoothforecast=Mfs),
class="STARS"
))
}, error = function(e) { stop(paste0("model fitting and forecasting are unsuccessful - please make sure the data and age range are suitable for the model and curve\n",e$message),call.=FALSE) })
}

#' @export
coef.STARS <- function(object,...) {
list(mu=object$mu,R=object$R)
}

#' @export
forecast.STARS <- function(object,...) {
object$smoothforecast
}

#' @export
plot.STARS <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(1,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
colband <- colorRampPalette(c("black","grey","white"))
image(x=x$x,y=c(2:nrow(x$M)),z=t(x$standardresiduals),col=colband(6),xlab="age",ylab="year",main="standardised residuals",cex.main=1,font.main=1)
plot(x$x,log(x$M[1,]),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n",ylim=c(min(log(x$M),log(x$smoothforecast)),max(log(x$M),log(x$smoothforecast))))
points(x$x,log(x$M[nrow(x$M),]),pch=1,cex=0.5)
lines(x$x,log(x$smoothforecast[x$h,]))
temp <- paste("forecast",x$h,"years")
legend("bottomright",legend=c("observed first data","observed last data",temp),pch=c(16,1,NA),lty=c(NA,NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
residuals.STARS <- function(object,...) {
object$standardresiduals
}
