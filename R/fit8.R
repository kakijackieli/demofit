#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit8 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x)))^2) }
suppressWarnings(resulta <- nlminb(c(1,0.1,0.001,1,120),f))
suppressWarnings(resultb <- optim(c(1,0.1,0.001,1,120),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x))) }
suppressWarnings(resultc <- nls.lm(c(1,0.1,0.001,1,120),h,lower=c(-Inf,-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf,Inf)))
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
I <- resulta$par[4]
N <- resulta$par[5]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
I <- resultb$par[4]
N <- resultb$par[5]
} else if (ind==3) {
A <- resultc$par[1]
B <- resultc$par[2]
C <- resultc$par[3]
I <- resultc$par[4]
N <- resultc$par[5]
}
fitted <- A+B*x+C*x^2+I/(N-x)
structure(
list(curve="Van der Maen",x=x,m=m,w=w,A=A,B=B,C=C,I=I,N=N,fitted=fitted,diagnostics=diagnostics),
class="Fit8"
)
}

#' @export
coef.Fit8 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,I=object$I,N=object$N)
}

#' @export
fitted.Fit8 <- .fittedFit

#' @export
predict.Fit8 <- function(object,newdata,...) {
object$A+object$B*newdata+object$C*newdata^2+object$I/(object$N-newdata)
}

#' @export
plot.Fit8 <- .plotFit

#' @export
deviance.Fit8 <- .devianceFit

#' @export
residuals.Fit8 <- .residualsFit
