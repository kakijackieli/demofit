#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef nlminb
#' @importFrom graphics lines legend

cfit7 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<=0)) { stop("x must be positive") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
logm <- log(m)
logx <- log(x)
fit <- lm(logm~logx,weights=w)
B <- as.numeric(exp(coef(fit)[1]))
C <- as.numeric(coef(fit)[2])
if (C<0) {
f <- function(p) { sum(w *(log(m)-log(p[1]*x^p[2]))^2) }
suppressWarnings(resulta <- nlminb(c(0.000000000001,10),f,lower=c(0,0),upper=c(Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]*x^p[2])) }
suppressWarnings(resultb <- nls.lm(c(0.000000000001,10),h,lower=c(0,0),upper=c(Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob))
if (ind==1) {
B <- resulta$par[1]
C <- resulta$par[2]
} else if (ind==2) {
B <- resultb$par[1]
C <- resultb$par[2]
}}
fitted <- B*x^C
structure(
list(curve="Weibull",x=x,m=m,w=w,B=B,C=C,fitted=fitted),
class="cFit7"
)
}

#' @export
coef.cFit7 <- function(object,...) {
c(B=object$B,C=object$C)
}

#' @export
fitted.cFit7 <- function(object,...) {
object$fitted
}

#' @export
predict.cFit7 <- function(object,newdata,...) {
object$B*newdata^object$C
}

#' @export
plot.cFit7 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.cFit7 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.cFit7 <- function(object,...) {
log(object$m)-log(object$fitted)
}
