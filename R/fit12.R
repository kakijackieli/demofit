#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit12 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]/exp(p[2]*x)+p[3]+p[4]*exp(p[5]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.01,10,0.0001,0.0001,0.1),f))
suppressWarnings(resultb <- optim(c(0.01,10,0.0001,0.0001,0.1),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]/exp(p[2]*x)+p[3]+p[4]*exp(p[5]*x))) }
suppressWarnings(resultc <- nls.lm(c(0.01,10,0.0001,0.0001,0.1),h,lower=c(-Inf,-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,nelder=resultb,levenberg=resultc)
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A1 <- resulta$par[1]
B1 <- resulta$par[2]
A2 <- resulta$par[3]
A3 <- resulta$par[4]
B3 <- resulta$par[5]
} else if (ind==2) {
A1 <- resultb$par[1]
B1 <- resultb$par[2]
A2 <- resultb$par[3]
A3 <- resultb$par[4]
B3 <- resultb$par[5]
} else if (ind==3) {
A1 <- resultc$par[1]
B1 <- resultc$par[2]
A2 <- resultc$par[3]
A3 <- resultc$par[4]
B3 <- resultc$par[5]
}
fitted <- A1/exp(B1*x)+A2+A3*exp(B3*x)
structure(
list(curve="Siler",x=x,m=m,w=w,A1=A1,B1=B1,A2=A2,A3=A3,B3=B3,fitted=fitted,diagnostics=diagnostics),
class="Fit12"
)
}

#' @export
coef.Fit12 <- function(object,...) {
c(A1=object$A1,B1=object$B1,A2=object$A2,A3=object$A3,B3=object$B3)
}

#' @export
fitted.Fit12 <- .fittedFit

#' @export
predict.Fit12 <- function(object,newdata,...) {
object$A1/exp(object$B1*newdata)+object$A2+object$A3*exp(object$B3*newdata)
}

#' @export
plot.Fit12 <- .plotFit

#' @export
deviance.Fit12 <- .devianceFit

#' @export
residuals.Fit12 <- .residualsFit
