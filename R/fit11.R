#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit11 <- function(x,m,w) {
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
f2 <- function(p) { sum(w[x>=10&x<=29]*(log(m[x>=10&x<=29])-log(p[1]/exp(p[2]*(x[x>=10&x<=29]-p[4])+1/exp(p[3]*(x[x>=10&x<=29]-p[4])))))^2) }
suppressWarnings(result2 <- nlminb(c(0.0001,0.1,0.1,20),f2))
f3 <- function(p) { sum(w[x>=30]*(log(m[x>=30])-log(p[1]*exp(p[2]*x[x>=30])))^2) }
suppressWarnings(result3 <- nlminb(c(0.0001,0.1),f3))
param = c(0.0001,result1$par[1],result2$par[1],result3$par[1],result1$par[2],result2$par[2:3],result3$par[2],result2$par[4])
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]/exp(p[5]*x)+p[3]/exp(p[6]*(x-p[9])+1/exp(p[7]*(x-p[9])))+p[4]*exp(p[8]*x)))^2) }
suppressWarnings(resulta <- nlminb(param,f))
suppressWarnings(resultb <- optim(param,f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]/exp(p[5]*x)+p[3]/exp(p[6]*(x-p[9])+1/exp(p[7]*(x-p[9])))+p[4]*exp(p[8]*x))) }
suppressWarnings(resultc <- nls.lm(param,h,lower=rep(-Inf,9),upper=rep(Inf,9)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(resultb$value),resultb$value,Inf)
oc = ifelse (is.finite(sum(resultc$fvec^2)),sum(resultc$fvec^2),Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
ind = which.min(c(oa,ob,oc))
if (ind==1) {
A0 <- resulta$par[1]
A1 <- resulta$par[2]
A2 <- resulta$par[3]
A3 <- resulta$par[4]
A <- resulta$par[5]
B <- resulta$par[6]
C <- resulta$par[7]
D <- resulta$par[8]
U <- resulta$par[9]
} else if (ind==2) {
A0 <- resultb$par[1]
A1 <- resultb$par[2]
A2 <- resultb$par[3]
A3 <- resultb$par[4]
A <- resultb$par[5]
B <- resultb$par[6]
C <- resultb$par[7]
D <- resultb$par[8]
U <- resultb$par[9]
} else if (ind==3) {
A0 <- resultc$par[1]
A1 <- resultc$par[2]
A2 <- resultc$par[3]
A3 <- resultc$par[4]
A <- resultc$par[5]
B <- resultc$par[6]
C <- resultc$par[7]
D <- resultc$par[8]
U <- resultc$par[9]
}
fitted <- A0+A1/exp(A*x)+A2/exp(B*(x-U)+1/exp(C*(x-U)))+A3*exp(D*x)
structure(
list(curve="Rogers Planck",x=x,m=m,w=w,A0=A0,A1=A1,A2=A2,A3=A3,A=A,B=B,C=C,D=D,U=U,fitted=fitted),
class="Fit11"
)
}

#' @export
coef.Fit11 <- function(object,...) {
c(A0=object$A0,A1=object$A1,A2=object$A2,A3=object$A3,A=object$A,B=object$B,C=object$C,D=object$D,U=object$U)
}

#' @export
fitted.Fit11 <- function(object,...) {
object$fitted
}

#' @export
predict.Fit11 <- function(object,newdata,...) {
object$A0+object$A1/exp(object$A*newdata)+object$A2/exp(object$B*(newdata-object$U)+1/exp(object$C*(newdata-object$U)))+object$A3*exp(object$D*newdata)
}

#' @export
plot.Fit11 <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

#' @export
deviance.Fit11 <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

#' @export
residuals.Fit11 <- function(object,...) {
log(object$m)-log(object$fitted)
}
