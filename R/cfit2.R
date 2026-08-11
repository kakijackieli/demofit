#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit2 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*exp(p[3]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,0.1),f,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*exp(p[3]*x))) }
suppressWarnings(resultb <- nls.lm(c(0.001,0.00001,0.1),h,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
ind = which.min(c(oa,ob))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
C <- resulta$par[3]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
}
fitted <- A+B*exp(C*x)
structure(
list(curve="Makeham",x=x,m=m,w=w,A=A,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="cFit2"
)
}

#' @export
coef.cFit2 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C)
}

#' @export
fitted.cFit2 <- .fittedFit

#' @export
predict.cFit2 <- function(object,newdata,...) {
object$A+object$B*exp(object$C*newdata)
}

#' @export
plot.cFit2 <- .plotFit

#' @export
deviance.cFit2 <- .devianceFit

#' @export
residuals.cFit2 <- .residualsFit
