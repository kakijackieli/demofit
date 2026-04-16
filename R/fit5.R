#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit5 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
f <- function(p) { sum(w *(log(m)-log(1/p[1]^((p[2]*x)^p[4])/p[2]+1/p[1]^((p[3]-x)^p[4])))^2) }
suppressWarnings(resulta <- nlminb(c(1,100,100,1),f))
suppressWarnings(resultb <- optim(c(1,100,100,1),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(1/p[1]^((p[2]*x)^p[4])/p[2]+1/p[1]^((p[3]-x)^p[4]))) }
suppressWarnings(resultc <- nls.lm(c(1,100,100,1),h,lower=c(-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
M <- resulta$par[3]
N <- resulta$par[4]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
M <- resultb$par[3]
N <- resultb$par[4]
} else if (ind==3) {
A <- resultc$par[1]
B <- resultc$par[2]
M <- resultc$par[3]
N <- resultc$par[4]
}
fitted <- 1/A^((B*x)^N)/B+1/A^((M-x)^N)
structure(
list(curve="Wittstein Bumsted",x=x,m=m,w=w,A=A,B=B,M=M,N=N,fitted=fitted),
class="Fit5"
)
}

#' @export
coef.Fit5 <- function(object,...) {
c(A=object$A,B=object$B,M=object$M,N=object$N)
}

#' @export
fitted.Fit5 <- function(object,...) {
object$fitted
}

#' @export
predict.Fit5 <- function(object,newdata,...) {
1/object$A^((object$B*newdata)^object$N)/object$B+1/object$A^((object$M-newdata)^object$N)
}

#' @export
plot.Fit5 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.Fit5 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.Fit5 <- function(object,...) {
log(object$m)-log(object$fitted)
}
