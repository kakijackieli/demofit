#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats nlminb coef
#' @importFrom graphics lines legend

cfit9 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]*exp(p[3]*x)/(1+p[2]*exp(p[3]*x))))^2) }
suppressWarnings(resulta <- nlminb(c(0.00001,0.00001,0.1),f,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]*exp(p[3]*x)/(1+p[2]*exp(p[3]*x)))) }
suppressWarnings(resultb <- nls.lm(c(0.00001,0.00001,0.1),h,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
suppressWarnings(resultc <- MortalityLaw(x=x,mx=m,law="beard"))
p1 <- as.numeric(coef(resultc)[1])/exp(as.numeric(coef(resultc)[2])*(min(x)-1))
p2 <- as.numeric(coef(resultc)[1])*as.numeric(coef(resultc)[3])/exp(as.numeric(coef(resultc)[2])*(min(x)-1))
p3 <- as.numeric(coef(resultc)[2])
error <- f(c(p1,p2,p3))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
oc = ifelse (is.finite(error),error,Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
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
A <- p1
B <- p2
C <- p3
}
fitted <- A*exp(C*x)/(1+B*exp(C*x))
structure(
list(curve="Beard",x=x,m=m,w=w,A=A,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="cFit9"
)
}

#' @export
coef.cFit9 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C)
}

#' @export
fitted.cFit9 <- .fittedFit

#' @export
predict.cFit9 <- function(object,newdata,...) {
object$A*exp(object$C*newdata)/(1+object$B*exp(object$C*newdata))
}

#' @export
plot.cFit9 <- .plotFit

#' @export
deviance.cFit9 <- .devianceFit

#' @export
residuals.cFit9 <- .residualsFit
