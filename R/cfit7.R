#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef nlminb
#' @importFrom graphics lines legend

cfit7 <- function(x,m,w) {
.checkinput(x,m,w)
if (any(x<=0)) { stop("x must be positive") }
logm <- log(m)
logx <- log(x)
fit <- lm(logm~logx,weights=w)
B <- as.numeric(exp(coef(fit)[1]))
C <- as.numeric(coef(fit)[2])
diagnostics <- NULL
if (C<0) {
f <- function(p) { sum(w *(log(m)-log(p[1]*x^p[2]))^2) }
suppressWarnings(resulta <- nlminb(c(0.000000000001,10),f,lower=c(0,0),upper=c(Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]*x^p[2])) }
suppressWarnings(resultb <- nls.lm(c(0.000000000001,10),h,lower=c(0,0),upper=c(Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
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
list(curve="Weibull",x=x,m=m,w=w,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="cFit7"
)
}

#' @export
coef.cFit7 <- function(object,...) {
c(B=object$B,C=object$C)
}

#' @export
fitted.cFit7 <- .fittedFit

#' @export
predict.cFit7 <- function(object,newdata,...) {
object$B*newdata^object$C
}

#' @export
plot.cFit7 <- .plotFit

#' @export
deviance.cFit7 <- .devianceFit

#' @export
residuals.cFit7 <- .residualsFit
