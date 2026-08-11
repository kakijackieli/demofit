#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit14 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*exp(p[3]*x)/(1+p[4]*exp(p[3]*x))))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.0001,0.1,0.0001),f,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*exp(p[3]*x)/(1+p[4]*exp(p[3]*x)))) }
suppressWarnings(resultb <- nls.lm(c(0.001,0.0001,0.1,0.0001),h,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
ind = which.min(c(oa,ob))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
C <- resulta$par[3]
D <- resulta$par[4]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
D <- resultb$par[4]
}
fitted <- A+B*exp(C*x)/(1+D*exp(C*x))
structure(
list(curve="Thatcher",x=x,m=m,w=w,A=A,B=B,C=C,D=D,fitted=fitted,diagnostics=diagnostics),
class="cFit14"
)
}

#' @export
coef.cFit14 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D)
}

#' @export
fitted.cFit14 <- .fittedFit

#' @export
predict.cFit14 <- function(object,newdata,...) {
object$A+object$B*exp(object$C*newdata)/(1+object$D*exp(object$C*newdata))
}

#' @export
plot.cFit14 <- .plotFit

#' @export
deviance.cFit14 <- .devianceFit

#' @export
residuals.cFit14 <- .residualsFit
