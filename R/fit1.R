#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef
#' @importFrom graphics lines legend

fit1 <- function(x,m,w) {
.checkinput(x,m,w)
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
fitted.Fit1 <- .fittedFit

#' @export
predict.Fit1 <- function(object,newdata,...) {
object$B*exp(object$C*newdata)
}

#' @export
plot.Fit1 <- .plotFit

#' @export
deviance.Fit1 <- .devianceFit

#' @export
residuals.Fit1 <- .residualsFit
