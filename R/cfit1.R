#' @importFrom minpack.lm nls.lm
#' @importFrom stats lm coef nlminb
#' @importFrom graphics lines legend

cfit1 <- function(x,m,w) {
.checkinput(x,m,w)
logm <- log(m)
fit <- lm(logm~x,weights=w)
B <- as.numeric(exp(coef(fit)[1]))
C <- as.numeric(coef(fit)[2])
diagnostics <- NULL
if (C<0) {
f <- function(p) { sum(w *(log(m)-log(p[1]*exp(p[2]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.00001,0.1),f,lower=c(0,0),upper=c(Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]*exp(p[2]*x))) }
suppressWarnings(resultb <- nls.lm(c(0.00001,0.1),h,lower=c(0,0),upper=c(Inf,Inf)))
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
fitted <- B*exp(C*x)
structure(
list(curve="Gompertz",x=x,m=m,w=w,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="cFit1"
)
}

#' @export
coef.cFit1 <- function(object,...) {
c(B=object$B,C=object$C)
}

#' @export
fitted.cFit1 <- .fittedFit

#' @export
predict.cFit1 <- function(object,newdata,...) {
object$B*exp(object$C*newdata)
}

#' @export
plot.cFit1 <- .plotFit

#' @export
deviance.cFit1 <- .devianceFit

#' @export
residuals.cFit1 <- .residualsFit
