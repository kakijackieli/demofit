#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef
#' @importFrom graphics lines legend

fit7 <- function(x,m,w) {
.checkinput(x,m,w)
if (any(x<=0)) { stop("x must be positive") }
logm <- log(m)
logx <- log(x)
fit <- lm(logm~logx,weights=w)
B <- as.numeric(exp(coef(fit)[1]))
C <- as.numeric(coef(fit)[2])
fitted <- B*x^C
structure(
list(curve="Weibull",x=x,m=m,w=w,B=B,C=C,fitted=fitted),
class="Fit7"
)
}

#' @export
coef.Fit7 <- function(object,...) {
c(B=object$B,C=object$C)
}

#' @export
fitted.Fit7 <- .fittedFit

#' @export
predict.Fit7 <- function(object,newdata,...) {
object$B*newdata^object$C
}

#' @export
plot.Fit7 <- .plotFit

#' @export
deviance.Fit7 <- .devianceFit

#' @export
residuals.Fit7 <- .residualsFit
