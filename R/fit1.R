#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef
#' @importFrom graphics lines legend

fit1 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
logm <- log(m)
fit <- lm(logm~x,weights=w)
B <- as.numeric(exp(coef(fit)[1]))
C <- as.numeric(coef(fit)[2])
fitted <- B*exp(C*x)
structure(
list(curve="Gompertz",x=x,m=m,w=w,B=B,C=C,fitted=fitted),
class="Fit1"
)
}

#' @export
coef.Fit1 <- function(object,...) {
c(B=object$B,C=object$C)
}

#' @export
fitted.Fit1 <- function(object,...) {
object$fitted
}

#' @export
predict.Fit1 <- function(object,newdata,...) {
object$B*exp(object$C*newdata)
}

#' @export
plot.Fit1 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.Fit1 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.Fit1 <- function(object,...) {
log(object$m)-log(object$fitted)
}
