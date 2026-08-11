#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats nlminb coef
#' @importFrom graphics lines legend

cfit10 <- function(x,m,w) {
.checkinput(x,m,w)
if ((min(x)>5)||(max(x)<60)) { stop("youngest age must be 5 or lower and oldest age must be 60 or higher") }
f1 <- function(p) { sum(w[x<=9]*(log(m[x<=9])-log(p[1]^((x[x<=9]+p[2])^p[3])))^2) }
suppressWarnings(result1 <- nlminb(c(0.001,0.01,0.1),f1,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
f2 <- function(p) { sum(w[x>=10&x<=29]*(log(m[x>=10&x<=29])-log(p[1]/exp(p[2]*(log(x[x>=10&x<=29])-log(p[3]))^2)))^2) }
suppressWarnings(result2 <- nlminb(c(0.001,10,20),f2,lower=c(0,0,0),upper=c(Inf,Inf,Inf)))
f3 <- function(p) { sum(w[x>=30]*(log(m[x>=30])-log(p[1]*p[2]^x[x>=30]/(1+p[1]*p[2]^x[x>=30])))^2) }
suppressWarnings(result3 <- nlminb(c(0.0001,1.1),f3,lower=c(0,0),upper=c(Inf,Inf)))
xx <- ifelse(x==0,1e-10,x)
f <- function(p) { sum(w *(log(m)-log(p[1]^((x+p[2])^p[3])+p[4]/exp(p[5]*(log(xx)-log(p[6]))^2)+p[7]*p[8]^x/(1+p[7]*p[8]^x)))^2) }
suppressWarnings(resulta <- nlminb(c(result1$par,result2$par,result3$par),f,lower=rep(0,8),upper=rep(Inf,8)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]^((x+p[2])^p[3])+p[4]/exp(p[5]*(log(xx)-log(p[6]))^2)+p[7]*p[8]^x/(1+p[7]*p[8]^x))) }
suppressWarnings(resultb <- nls.lm(c(result1$par,result2$par,result3$par),h,lower=rep(0,8),upper=rep(Inf,8)))
resultc <- MortalityLaw(x=x,mx=m,law="HP2")
error <- f(as.numeric(coef(resultc)))
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
D <- resulta$par[4]
E <- resulta$par[5]
F <- resulta$par[6]
G <- resulta$par[7]
H <- resulta$par[8]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
D <- resultb$par[4]
E <- resultb$par[5]
F <- resultb$par[6]
G <- resultb$par[7]
H <- resultb$par[8]
} else if (ind==3) {
A <- as.numeric(coef(resultc)[1])
B <- as.numeric(coef(resultc)[2])
C <- as.numeric(coef(resultc)[3])
D <- as.numeric(coef(resultc)[4])
E <- as.numeric(coef(resultc)[5])
F <- as.numeric(coef(resultc)[6])
G <- as.numeric(coef(resultc)[7])
H <- as.numeric(coef(resultc)[8])
}
fitted <- A^((x+B)^C)+D/exp(E*(log(xx)-log(F))^2)+G*H^x/(1+G*H^x)
structure(
list(curve="Heligman Pollard",x=x,m=m,w=w,A=A,B=B,C=C,D=D,E=E,F=F,G=G,H=H,fitted=fitted,diagnostics=diagnostics),
class="cFit10"
)
}

#' @export
coef.cFit10 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D,E=object$E,F=object$F,G=object$G,H=object$H)
}

#' @export
fitted.cFit10 <- .fittedFit

#' @export
predict.cFit10 <- function(object,newdata,...) {
xx <- ifelse(newdata==0,1e-10,newdata)
object$A^((newdata+object$B)^object$C)+object$D/exp(object$E*(log(xx)-log(object$F))^2)+object$G*object$H^newdata/(1+object$G*object$H^newdata)
}

#' @export
plot.cFit10 <- .plotFit

#' @export
deviance.cFit10 <- .devianceFit

#' @export
residuals.cFit10 <- .residualsFit
