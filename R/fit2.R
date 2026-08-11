#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb 
#' @importFrom graphics lines legend

fit2 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*exp(p[3]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,0.1),f))
suppressWarnings(resultb <- optim(c(0.001,0.00001,0.1),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*exp(p[3]*x))) }
suppressWarnings(resultc <- nls.lm(c(0.001,0.00001,0.1),h,lower=c(-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,nelder=resultb,levenberg=resultc)
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
C <- resulta$par[3]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
} else if (ind==3) {
A <- resultc$par[1]
B <- resultc$par[2]
C <- resultc$par[3]
}
fitted <- A+B*exp(C*x)
structure(
list(curve="Makeham",x=x,m=m,w=w,A=A,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="Fit2"
)
}

#' @export
coef.Fit2 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C)
}

#' @export
fitted.Fit2 <- .fittedFit

#' @export
predict.Fit2 <- function(object,newdata,...) {
object$A+object$B*exp(object$C*newdata)
}

#' @export
plot.Fit2 <- .plotFit

#' @export
deviance.Fit2 <- .devianceFit

#' @export
residuals.Fit2 <- .residualsFit
