#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit3 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]/sqrt(x+1)-p[2]+p[3]*sqrt(x+1)))^2) }
suppressWarnings(resulta <- nlminb(c(0.01,0.001,0.001),f,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]/sqrt(x+1)-p[2]+p[3]*sqrt(x+1))) }
suppressWarnings(resultb <- nls.lm(c(0.01,0.001,0.001),h,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
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
fitted <- A/sqrt(x+1)-B+C*sqrt(x+1)
structure(
list(curve="Oppermann",x=x,m=m,w=w,A=A,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="cFit3"
)
}

#' @export
coef.cFit3 <- function(object,...) {
c(A=object$A,B=-object$B,C=object$C)
}

#' @export
fitted.cFit3 <- .fittedFit

#' @export
predict.cFit3 <- function(object,newdata,...) {
object$A/sqrt(newdata+1)-object$B+object$C*sqrt(newdata+1)
}

#' @export
plot.cFit3 <- .plotFit

#' @export
deviance.cFit3 <- .devianceFit

#' @export
residuals.cFit3 <- .residualsFit
