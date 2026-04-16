#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit12 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
f <- function(p) { sum(w *(log(m)-log(p[1]/exp(p[2]*x)+p[3]+p[4]*exp(p[5]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.01,10,0.0001,0.0001,0.1),f,lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]/exp(p[2]*x)+p[3]+p[4]*exp(p[5]*x))) }
suppressWarnings(resultb <- nls.lm(c(0.01,10,0.0001,0.0001,0.1),h,lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob))
if (ind==1) {
A1 <- resulta$par[1]
B1 <- resulta$par[2]
A2 <- resulta$par[3]
A3 <- resulta$par[4]
B3 <- resulta$par[5]
} else if (ind==2) {
A1 <- resultb$par[1]
B1 <- resultb$par[2]
A2 <- resultb$par[3]
A3 <- resultb$par[4]
B3 <- resultb$par[5]
}
fitted <- A1/exp(B1*x)+A2+A3*exp(B3*x)
structure(
list(curve="Siler",x=x,m=m,w=w,A1=A1,B1=B1,A2=A2,A3=A3,B3=B3,fitted=fitted),
class="cFit12"
)
}

#' @export
coef.cFit12 <- function(object,...) {
c(A1=object$A1,B1=object$B1,A2=object$A2,A3=object$A3,B3=object$B3)
}

#' @export
fitted.cFit12 <- function(object,...) {
object$fitted
}

#' @export
predict.cFit12 <- function(object,newdata,...) {
object$A1/exp(object$B1*newdata)+object$A2+object$A3*exp(object$B3*newdata)
}

#' @export
plot.cFit12 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.cFit12 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.cFit12 <- function(object,...) {
log(object$m)-log(object$fitted)
}
