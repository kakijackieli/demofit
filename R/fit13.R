#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats optim nlminb coef
#' @importFrom graphics lines legend

fit13 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
f <- function(p) { sum(w *(log(m)-log((p[1]+p[2]*exp(p[3]*x))/(1+p[4]*exp(p[3]*x))+p[5]*exp(p[3]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,0.1,0.00001,0.0001),f))
suppressWarnings(resultb <- optim(c(0.001,0.00001,0.1,0.00001,0.0001),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log((p[1]+p[2]*exp(p[3]*x))/(1+p[4]*exp(p[3]*x))+p[5]*exp(p[3]*x))) }
suppressWarnings(resultc <- nls.lm(c(0.001,0.00001,0.1,0.00001,0.0001),h,lower=c(-Inf,-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf,Inf)))
resultd <- MortalityLaw(x=x,mx=m,law="martinelle")
error <- f(as.numeric(c(coef(resultd)[3],coef(resultd)[1],coef(resultd)[2],coef(resultd)[4],coef(resultd)[5])))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
od = ifelse (is.finite(error),error,Inf)
if (all(!is.finite(c(oa,ob,oc,od)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob,oc,od))
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
A <- resultc$par[1]
B <- resultc$par[2]
C <- resultc$par[3]
D <- resultc$par[4]
E <- resultc$par[5]
} else if (ind==4) {
A <- as.numeric(coef(resultd)[3])
B <- as.numeric(coef(resultd)[1])
C <- as.numeric(coef(resultd)[2])
D <- as.numeric(coef(resultd)[4])
E <- as.numeric(coef(resultd)[5])
}
fitted <- (A+B*exp(C*x))/(1+D*exp(C*x))+E*exp(C*x)
structure(
list(curve="Martinelle",x=x,m=m,w=w,A=A,B=B,C=C,D=D,E=E,fitted=fitted),
class="Fit13"
)
}

#' @export
coef.Fit13 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D,E=object$E)
}

#' @export
fitted.Fit13 <- function(object,...) {
object$fitted
}

#' @export
predict.Fit13 <- function(object,newdata,...) {
(object$A+object$B*exp(object$C*newdata))/(1+object$D*exp(object$C*newdata))+object$E*exp(object$C*newdata)
}

#' @export
plot.Fit13 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.Fit13 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.Fit13 <- function(object,...) {
log(object$m)-log(object$fitted)
}
