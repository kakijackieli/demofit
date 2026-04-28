#' Mortality ensemble interval
#'
#' Generates ensemble interval forecast of mortality rates.
#'
#' @param ... fitted model objects returned by \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, and / or \code{STARS()}.
#' @param width coverage probability of interval (default = 0.95).
#' @param method if 1, simple averaging; if 2, weighted averaging; if 3, envelope; if 4, interior trimming (default = 1).
#' @param wm vector of weights for LC, RH, APC, M5, M6, M7, and / or STAR models if method = 2 (default = equal weights).
#' @param nsim number of simulations (default = 10).
#' @param seed seed for random number generator (default = 123).
#'
#' @details
#' Ensemble interval forecast is constructed by combining interval forecasts from individual stochastic mortality models using different methods including simple averaging, weighted averaging, envelope, and interior trimming. See \code{LCS()}, \code{RHS()}, \code{APCS()}, \code{CBDS()}, \code{CBDCS()}, \code{CBDQCS()}, and \code{STARS()} for more details of different stochastic mortality models.
#'
#' @importFrom stats simulate quantile
#'
#' @return 
#' An object of class ENI with associated S3 methods forecast and plot.
#'
#' @references
#' Li, J., Wang, M., Liu, J., and Tickle, L. (2025). Ensemble interval forecasts of mortality. Scandinavian Actuarial Journal, 2025(6), 598-616.
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
#' fit1 <- LCS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' fit2 <- RHS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' fit3 <- APCS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' fit4 <- CBDS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
#' fit <- ENI(fit1,fit2,fit3,fit4)
#' forecast::forecast(fit)
#' plot(fit)
#'
#' @export
ENI <- function(...,width=0.95,method=1,wm,nsim=10,seed=123) {
fits <- list(...)
n <- length(fits)
if (n<2||n>7) { stop("number of models must be between 2 and 7") }
if (missing(wm)) { wm <- rep(1/n,n) }
allowed <- c("LCS","APCS","RHS","CBDS","CBDCS","CBDQCS","STARS")
classes <- sapply(fits,function(x) class(x)[1])
if (!all(classes %in% allowed)) { stop("all inputs must be valid fitted model objects") }
if (any(duplicated(classes))) { stop("fitted model objects must be distinct, i.e. no duplicates") }
x <- fits[[1]]$x
h <- fits[[1]]$h  
for (f in fits) {
if (!identical(f$x,x)) stop("all models must have same x")
if (f$h!=h) stop("all models must have same h") }
if (!is.numeric(width)||width<=0||width>=1) { stop("width must be a number between 0 and 1") }
if (!is.numeric(method)||!(method %in% 1:4)) { stop("method must be 1, 2, 3, or 4") }
if (method==2) {
if (length(wm)!=n) { stop("length of wm must be equal to number of models") }
if (any(wm<0)) { stop("wm must be non-negative") }
if (abs(sum(wm)-1)>1e-8) { stop("wm must sum to 1") }}
tryCatch({
nc <- length(x)
sims <- lapply(fits,function(f) simulate(f,nsim=nsim,seed=seed))
if (method==1) {
upper <- matrix(0,h,nc); lower <- matrix(0,h,nc); z <- 1
for (f in fits) {
upper <- upper+apply(sims[[z]]$Msim,c(1,2),quantile,probs=1-(1-width)/2)/n
lower <- lower+apply(sims[[z]]$Msim,c(1,2),quantile,probs=(1-width)/2)/n 
z <- z+1 }}
if (method==2) { 
upper <- matrix(0,h,nc); lower <- matrix(0,h,nc); z <- 1
for (f in fits) {
upper <- upper+apply(sims[[z]]$Msim,c(1,2),quantile,probs=1-(1-width)/2)*wm[z]
lower <- lower+apply(sims[[z]]$Msim,c(1,2),quantile,probs=(1-width)/2)*wm[z] 
z <- z+1 }}
if (method==3) {
upper <- matrix(-Inf,h,nc); lower <- matrix(Inf,h,nc); z <- 1
for (f in fits) {
upper <- pmax(upper,apply(sims[[z]]$Msim,c(1,2),quantile,probs=1-(1-width)/2))
lower <- pmin(lower,apply(sims[[z]]$Msim,c(1,2),quantile,probs=(1-width)/2)) 
z <- z+1 }}
if (method==4) {
u <- matrix(Inf,h,nc); l <- matrix(-Inf,h,nc); z <- 1
for (f in fits) {
u <- pmin(u,apply(sims[[z]]$Msim,c(1,2),quantile,probs=1-(1-width)/2))
l <- pmax(l,apply(sims[[z]]$Msim,c(1,2),quantile,probs=(1-width)/2)) 
z <- z+1 }
upper <- matrix(0,h,nc); lower <- matrix(0,h,nc); z <- 1
for (f in fits) {
upper <- upper+apply(sims[[z]]$Msim,c(1,2),quantile,probs=1-(1-width)/2)/(n-1)
lower <- lower+apply(sims[[z]]$Msim,c(1,2),quantile,probs=(1-width)/2)/(n-1) 
z <- z+1 }
upper <- upper-u/(n-1)
lower <- lower-l/(n-1) }
invisible(structure(
list(x=x,M=fits[[1]]$M,wm=wm,h=h,upper=upper,lower=lower),
class="ENI"
))
}, error = function(e) { stop(paste0("ensemble is unsuccessful - please make sure the fitted model objects are valid\n",e$message),call.=FALSE) })
}

#' @export
forecast.ENI <- function(object,...) {
list(upper=object$upper,lower=object$lower)
}

#' @export
plot.ENI <- function(x,...) {
old <- par(no.readonly=TRUE)
on.exit(par(old))
par(mfrow=c(1,2),mar=c(3.5,2.5,1.5,0.5),mgp=c(1.5,0.5,0))
nr <- nrow(x$M); nc <- ncol(x$M)
plot(1:nr,log(x$M[,1]),xlab="year",ylab="log death rate (first age)",type="l",bty="n",xlim=c(1,nr+x$h),ylim=c(min(log(x$M[,1]),log(x$lower[,1])),max(log(x$M[,1]),log(x$upper[,1]))))
lines(nr+(1:x$h),log(x$upper[,1]),lty=2)
lines(nr+(1:x$h),log(x$lower[,1]),lty=2)
legend("bottomright",legend=c("observed","ensemble interval"),lty=c(1,2),cex=0.8,bty="n")
plot(1:nr,log(x$M[,nc]),xlab="year",ylab="log death rate (last age)",type="l",bty="n",xlim=c(1,nr+x$h),ylim=c(min(log(x$M[,nc]),log(x$lower[,nc])),max(log(x$M[,nc]),log(x$upper[,nc]))))
lines(nr+(1:x$h),log(x$upper[,nc]),lty=2)
lines(nr+(1:x$h),log(x$lower[,nc]),lty=2)
legend("bottomright",legend=c("observed","ensemble interval"),lty=c(1,2),cex=0.8,bty="n")
}
