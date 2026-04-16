#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit8 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x)))^2) }
suppressWarnings(resulta <- nlminb(c(1,0.1,0.001,1,120),f,lower=c(0,-Inf,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x))) }
suppressWarnings(resultb <- nls.lm(c(1,0.1,0.001,1,120),h,lower=c(0,-Inf,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
C <- resulta$par[3]
I <- resulta$par[4]
N <- resulta$par[5]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
I <- resultb$par[4]
N <- resultb$par[5]
}
fitted <- A+B*x+C*x^2+I/(N-x)
structure(
list(curve="Van der Maen",x=x,m=m,w=w,A=A,B=B,C=C,I=I,N=N,fitted=fitted),
class="cFit8"
)
}

#' @export
coef.cFit8 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,I=object$I,N=object$N)
}

#' @export
fitted.cFit8 <- function(object,...) {
object$fitted
}

#' @export
predict.cFit8 <- function(object,newdata,...) {
object$A+object$B*newdata+object$C*newdata^2+object$I/(object$N-newdata)
}

#' @export
plot.cFit8 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.cFit8 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.cFit8 <- function(object,...) {
log(object$m)-log(object$fitted)
}
