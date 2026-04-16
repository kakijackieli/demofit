#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit4 <- function(x,m,w) {
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
if ((min(x)>5)||(max(x)<60)) { stop("youngest age must be 5 or lower and oldest age must be 60 or higher") }
f1 <- function(p) { sum(w[x<=9]*(log(m[x<=9])-log(p[1]/exp(p[2]*x[x<=9])))^2) }
suppressWarnings(result1 <- nlminb(c(0.01,1),f1))
f2 <- function(p) { sum(w[x>=10&x<=29]*(log(m[x>=10&x<=29])-log(p[1]/exp(0.5*p[2]*(x[x>=10&x<=29]-p[3])^2)))^2) }
suppressWarnings(result2 <- nlminb(c(0.001,0.1,20),f2))
f3 <- function(p) { sum(w[x>=30]*(log(m[x>=30])-log(p[1]*exp(p[2]*x[x>=30])))^2) }
suppressWarnings(result3 <- nlminb(c(0.0001,0.1),f3))
f <- function(p) { sum(w *(log(m)-log(p[1]/exp(p[2]*x)+p[3]/exp(0.5*p[4]*(x-p[5])^2)+p[6]*exp(p[7]*x)))^2) }
suppressWarnings(resulta <- nlminb(c(result1$par,result2$par,result3$par),f))
suppressWarnings(resultb <- optim(c(result1$par,result2$par,result3$par),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]/exp(p[2]*x)+p[3]/exp(0.5*p[4]*(x-p[5])^2)+p[6]*exp(p[7]*x))) }
suppressWarnings(resultc <- nls.lm(c(result1$par,result2$par,result3$par),h,lower=rep(-Inf,7),upper=rep(Inf,7)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A1 <- resulta$par[1]
B1 <- resulta$par[2]
A2 <- resulta$par[3]
B2 <- resulta$par[4]
C <- resulta$par[5]
A3 <- resulta$par[6]
B3 <- resulta$par[7]
} else if (ind==2) {
A1 <- resultb$par[1]
B1 <- resultb$par[2]
A2 <- resultb$par[3]
B2 <- resultb$par[4]
C <- resultb$par[5]
A3 <- resultb$par[6]
B3 <- resultb$par[7]
} else if (ind==3) {
A1 <- resultc$par[1]
B1 <- resultc$par[2]
A2 <- resultc$par[3]
B2 <- resultc$par[4]
C <- resultc$par[5]
A3 <- resultc$par[6]
B3 <- resultc$par[7]
}
fitted <- A1/exp(B1*x)+A2/exp(0.5*B2*(x-C)^2)+A3*exp(B3*x)
structure(
list(curve="Thiele",x=x,m=m,w=w,A1=A1,B1=B1,A2=A2,B2=B2,C=C,A3=A3,B3=B3,fitted=fitted),
class="Fit4"
)
}

#' @export
coef.Fit4 <- function(object,...) {
c(A1=object$A1,B1=object$B1,A2=object$A2,B2=object$B2,C=object$C,A3=object$A3,B3=object$B3)
}

#' @export
fitted.Fit4 <- function(object,...) {
object$fitted
}

#' @export
predict.Fit4 <- function(object,newdata,...) {
object$A1/exp(object$B1*newdata)+object$A2/exp(0.5*object$B2*(newdata-object$C)^2)+object$A3*exp(object$B3*newdata)
}

#' @export
plot.Fit4 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.Fit4 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.Fit4 <- function(object,...) {
log(object$m)-log(object$fitted)
}
