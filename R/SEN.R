#' Lee-Carter model sensitivities
#'
#' Calculates sensitivities of Lee-Carter model parameters, mortality forecasts, and derived actuarial quantities to small changes in data.
#'
#' @param M matrix of central death rates (rows as years, columns as ages).
#' @param x index age of data point (column number of M) being perturbed.
#' @param t index year of data point (row number of M) being perturbed.
#' @param x1 index age of output variable (or initial age for annuity or assurance) whose sensitivity is calculated.
#' @param t1 index year of output variable whose sensitivity is calculated.
#' @param n1 term of annuity or assurance whose sensitivity is calculated.
#' @param r annual interest rate for annuity or assurance whose sensitivity is calculated.
#' @param f type of output variable whose sensitivity is calculated (including s, u, v, a, b, k, kf, mf, axn, A1xn).
#' @param jumpoff if 1, forecasts are based on estimated parameters only; if 2, forecasts are anchored to observed death rates in final year (default = 1). 
#' @param mq if 1, sensitivity is calculated with respect to central death rate; if 2, sensitivity is calculated with respect to probability of death (default = 1). 
#'
#' @details
#' User can choose one of the 10 output variables, including the first singular value (s), first left singular vector (u), first right singular vector (v), alpha (a), beta (b), kappa (k), forecasted kappa (kf), forecasted death rates (mf), expected present value of a term annuity ($1 payable at the end of each year of survival) (axn), and expected present value of a term assurance ($1 payable at the end of the year of death) (A1xn), under the Lee-Carter model, forecasted by applying the random walk with drift to \eqn{\kappa_t}. Their sensitivities are analytically derived using first-order matrix perturbation theory. 
#'
#' @importFrom forecast forecast
#'
#' @return
#' An object of class SEN with associated S3 methods coef and forecast.
#'
#' @references
#' Li, J. and Haberman, S. (2026). The Lee-Carter model revisited: Perturbation analysis and estimation equivalence. Working Paper.
#'
#' @examples
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
#' SEN(M,x=10,t=15,f="s")$sensitivity
#' SEN(M,x=10,t=15,x1=10,f="u")$sensitivity
#' SEN(M,x=10,t=15,t1=15,f="v")$sensitivity
#' SEN(M,x=10,t=15,x1=10,f="a")$sensitivity
#' SEN(M,x=10,t=15,x1=10,f="b")$sensitivity
#' SEN(M,x=10,t=15,t1=15,f="k")$sensitivity
#' SEN(M,x=10,t=15,t1=35,f="kf")$sensitivity
#' SEN(M,x=10,t=15,x1=10,t1=35,f="mf")$sensitivity
#' SEN(M,x=10,t=15,x1=10,n1=10,r=0.02,f="axn")$sensitivity
#' SEN(M,x=10,t=15,x1=10,n1=10,r=0.02,f="A1xn")$sensitivity
#'
#' @export
SEN <- function(M,x,t,x1=NULL,t1=NULL,n1=NULL,r=NULL,f=c("s","u","v","a","b","k","kf","mf","axn","A1xn"),jumpoff=1,mq=1) {
if (!is.numeric(M)) stop("M must be numeric")
if (!is.matrix(M)) stop("M must be a matrix with its rows as years and columns as ages")
if (any(!is.finite(M))) stop("all M values must be finite")
if (any(M<=0)) stop("all M values must be positive")
if (ncol(M)<2) stop("M must have at least 2 columns")
if (nrow(M)<2) stop("M must have at least 2 rows")
f <- tryCatch(match.arg(f),error = function(e) { stop("invalid variable choice") })
if (f %in% c("kf","mf","axn","A1xn")&&nrow(M)<20) stop("it requires at least 20 years of data for forecasting")
check <- function(z,name) { if (is.null(z)||!is.numeric(z)||length(z)!=1||!is.finite(z)||z<1||z!=floor(z)) stop(paste(name,"must be a positive integer")) }
check(x,"x")
check(t,"t")
if (f %in% c("u","a","b","mf","axn","A1xn")) check(x1,"x1")
if (f %in% c("v","k","kf","mf")) check(t1,"t1")
if (f %in% c("axn","A1xn")) {
check(n1,"n1")
if (is.null(r)||!is.numeric(r)||length(r)!=1||!is.finite(r)||r<0) stop("r must be a non-negative finite number")
}
if (!is.numeric(jumpoff)||length(jumpoff)!=1||!is.finite(jumpoff)||!(jumpoff %in% c(1,2))) stop("jumpoff must be either 1 or 2")
if (!is.numeric(mq)||length(mq)!=1||!is.finite(mq)||!(mq %in% c(1,2))) stop("mq must be either 1 or 2")
if (x>ncol(M)) stop("x must not exceed the number of columns of M")
if (t>nrow(M)) stop("t must not exceed the number of rows of M")
if (!is.null(x1)&&x1>ncol(M)) stop("x1 must not exceed the number of columns of M")
if (f %in% c("v","k")&&!is.null(t1)&&t1>nrow(M)) stop("t1 must not exceed the number of rows of M")
if (f %in% c("kf","mf")&&!is.null(t1)&&t1<=nrow(M)) stop("t1 must exceed the number of rows of M")
if (f %in% c("axn","A1xn")&&x1+n1-1>ncol(M)) stop("x1+n1-1 must not exceed the number of columns of M")
tryCatch({
Y <- t(log(M))
p <- nrow(Y); n <- ncol(Y)
a <- rowMeans(Y)
Z <- Y-a
S <- svd(Z,nu=p,nv=n)
s <- S$d
u <- S$u
v <- S$v
if (s[1]<=s[2]) stop("the leading singular value must be distinct")
if (abs(sum(u[,1]))<1e-10) stop("the first left singular vector must not sum to 0")
b <- u[,1]/sum(u[,1])
k <- s[1]*v[,1]*sum(u[,1])
M <- t(M)
kf <- k
Mf <- NULL
if (f %in% c("kf","mf","axn","A1xn")) {
if (f %in% c("kf","mf")) { h <- t1 }
if (f %in% c("axn","A1xn")) { h <- n+n1 }
Mf <- cbind(M,matrix(NA,nrow=p,ncol=h-n))
for (tt in (n+1):h) {
kf[tt] <- kf[tt-1]+(k[n]-k[1])/(n-1)
if (jumpoff==1) {
for (xx in 1:p) { Mf[xx,tt] <- exp(a[xx]+b[xx]*kf[tt]) }
} else {
for (xx in 1:p) { Mf[xx,tt] <- M[xx,n]*exp(b[xx]*(kf[tt]-k[n])) }
}}}
an <- function(x1,n1,r) {
temp1 <- 1; temp2 <- 0
for (tt in (n+1):(n+n1)) {
temp1 <- temp1/exp(Mf[x1+tt-n-1,tt])
temp2 <- temp2+temp1*(1+r)^(-tt+n)
}
temp2
}
as <- function(x1,n1,r) {
temp1 <- 1; temp2 <- 0
for (tt in (n+1):(n+n1)) {
temp2 <- temp2+(1-exp(-Mf[x1+tt-n-1,tt]))*temp1*(1+r)^(-tt+n)
temp1 <- temp1/exp(Mf[x1+tt-n-1,tt])
}
temp2
}
ds1 <- function(x,t) { u[x,1]*v[t,1]/M[x,t] }
du1 <- function(x1,x,t) {
temp1 <- 0; temp2 <- 0
for (i in 2:min(p,n)) { temp1 <- temp1+(s[1]*u[x,i]*v[t,1]+s[i]*u[x,1]*v[t,i])/(s[1]^2-s[i]^2)*u[x1,i] }
if (p>n) { for (i in (n+1):p) { temp2 <- temp2+u[x,i]*v[t,1]/s[1]*u[x1,i] }}
(temp1+temp2)/M[x,t]
}
dv1 <- function(t1,x,t) {
temp1 <- 0; temp2 <- 0; temp3 <- 0; temp4 <- 0
for (i in 2:min(p,n)) { temp1 <- temp1+(s[i]*u[x,i]*v[t,1]+s[1]*u[x,1]*v[t,i])/(s[1]^2-s[i]^2)*v[t1,i] }
if (n>p) { for (i in (p+1):n) { temp2 <- temp2+u[x,1]*v[t,i]/s[1]*v[t1,i] }}
for (ss in 1:n) {
for (i in 2:min(p,n)) { temp3 <- temp3+s[1]*u[x,1]*v[ss,i]/(s[1]^2-s[i]^2)*v[t1,i] }
if (n>p) { for (i in (p+1):n) { temp4 <- temp4+u[x,1]*v[ss,i]/s[1]*v[t1,i] }}
}
(temp1+temp2-temp3/n-temp4/n)/M[x,t]
}
da <- function(x1,x,t) {
if (x==x1) {
temp <- 1/n/M[x,t]
} else {
temp <- 0
}
temp
}
db <- function(x1,x,t) {
temp1 <- 0; temp2 <- 0
for (aa in 1:p) {
temp1 <- temp1+u[aa,1]*du1(x1,x,t)-u[x1,1]*du1(aa,x,t)
temp2 <- temp2+u[aa,1]
}
temp1/temp2^2
}
dk <- function(t1,x,t) {
temp1 <- 0; temp2 <- 0
for (aa in 1:p) {
temp1 <- temp1+u[aa,1]
temp2 <- temp2+du1(aa,x,t)
}
(ds1(x,t)*v[t1,1]+s[1]*dv1(t1,x,t))*temp1+s[1]*v[t1,1]*temp2
}
dkf <- function(t1,x,t) { (t1-1)/(n-1)*dk(n,x,t)-(t1-n)/(n-1)*dk(1,x,t) }
dmf <- function(x1,t1,x,t) {
if (jumpoff==1) {
temp <- Mf[x1,t1]*(da(x1,x,t)+kf[t1]*db(x1,x,t)+b[x1]*dkf(t1,x,t))
} else {
temp <- Mf[x1,t1]*((kf[t1]-k[n])*db(x1,x,t)+b[x1]*(dkf(t1,x,t)-dk(n,x,t)))
if (x==x1&&t==n) { temp <- temp+exp(b[x1]*(kf[t1]-k[n])) }
}
temp
}
dan <- function(x1,n1,x,t,r) {
temp1 <- 1; temp2 <- 0; temp3 <- 0
for (tt in (n+1):(n+n1)) {
temp1 <- temp1/exp(Mf[x1+tt-n-1,tt])
temp2 <- temp2+dmf(x1+tt-n-1,tt,x,t)
temp3 <- temp3-temp2*temp1*(1+r)^(-tt+n)
}
temp3
}
das <- function(x1,n1,x,t,r) {
temp1 <- 1; temp2 <- 0; temp3 <- 0
for (tt in (n+1):(n+n1)) {
temp3 <- temp3+(exp(-Mf[x1+tt-n-1,tt])*dmf(x1+tt-n-1,tt,x,t)*temp1-(1-exp(-Mf[x1+tt-n-1,tt]))*temp2*temp1)*(1+r)^(-tt+n)
temp1 <- temp1/exp(Mf[x1+tt-n-1,tt])
temp2 <- temp2+dmf(x1+tt-n-1,tt,x,t)
}
temp3
}
if (f=="s") sen <- ds1(x,t)
if (f=="u") sen <- du1(x1,x,t)
if (f=="v") sen <- dv1(t1,x,t)
if (f=="a") sen <- da(x1,x,t)
if (f=="b") sen <- db(x1,x,t)
if (f=="k") sen <- dk(t1,x,t)
if (f=="kf") sen <- dkf(t1,x,t)
if (f=="mf") sen <- dmf(x1,t1,x,t)
if (f=="axn") sen <- c(sensitivity=dan(x1,n1,x,t,r),EPV=an(x1,n1,r))
if (f=="A1xn") sen <- c(sensitivity=das(x1,n1,x,t,r),EPV=as(x1,n1,r))
if (mq==2) sen[1] <- sen[1]*exp(M[x,t])
invisible(structure(
list(M=M,x=x,t=t,x1=x1,t1=t1,n1=n1,r=r,f=f,jumpoff=jumpoff,mq=mq,sensitivity=sen,alpha=a,beta=b,kappa=k,forecast=Mf),
class="SEN"
))
}, error = function(e) { stop(paste0("sensitivity calculation is unsuccessful - please make sure the data structure is valid\n",e$message),call.=FALSE) })
}

#' @export
coef.SEN <- function(object,...) {
list(alpha=object$alpha,beta=object$beta,kappa=object$kappa)
}

#' @export
forecast.SEN <- function(object,...) {
object$forecast
}
