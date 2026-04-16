#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats nlminb coef
#' @importFrom graphics lines legend

cfit13 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
f <- function(p) { sum(w *(log(m)-log((p[1]+p[2]*exp(p[3]*x))/(1+p[4]*exp(p[3]*x))+p[5]*exp(p[3]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,0.1,0.00001,0.0001),f,lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log((p[1]+p[2]*exp(p[3]*x))/(1+p[4]*exp(p[3]*x))+p[5]*exp(p[3]*x))) }
suppressWarnings(resultb <- nls.lm(c(0.001,0.00001,0.1,0.00001,0.0001),h,lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
resultc <- MortalityLaw(x=x,mx=m,law="martinelle")
error <- f(as.numeric(c(coef(resultc)[3],coef(resultc)[1],coef(resultc)[2],coef(resultc)[4],coef(resultc)[5])))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
oc = ifelse (is.finite(error),error,Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
C <- resulta$par[3]
D <- resulta$par[4]
E <- resulta$par[5]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
D <- resultb$par[4]
E <- resultb$par[5]
} else if (ind==3) {
A <- as.numeric(coef(resultc)[3])
B <- as.numeric(coef(resultc)[1])
C <- as.numeric(coef(resultc)[2])
D <- as.numeric(coef(resultc)[4])
E <- as.numeric(coef(resultc)[5])
}
fitted <- (A+B*exp(C*x))/(1+D*exp(C*x))+E*exp(C*x)
structure(
list(curve="Martinelle",x=x,m=m,w=w,A=A,B=B,C=C,D=D,E=E,fitted=fitted),
class="cFit13"
)
}

#' @export
coef.cFit13 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D,E=object$E)
}

#' @export
fitted.cFit13 <- function(object,...) {
object$fitted
}

#' @export
predict.cFit13 <- function(object,newdata,...) {
(object$A+object$B*exp(object$C*newdata))/(1+object$D*exp(object$C*newdata))+object$E*exp(object$C*newdata)
}

#' @export
plot.cFit13 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.cFit13 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.cFit13 <- function(object,...) {
log(object$m)-log(object$fitted)
}
