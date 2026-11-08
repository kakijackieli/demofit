#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats nlminb coef
#' @importFrom graphics lines legend

cfit8 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x)))^2) }
suppressWarnings(resulta <- nlminb(c(1,0.1,0.001,1,120),f,lower=c(0,-Inf,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]+p[2]*x+p[3]*x^2+p[4]/(p[5]-x))) }
suppressWarnings(resultb <- nls.lm(c(1,0.1,0.001,1,120),h,lower=c(0,-Inf,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf)))
suppressWarnings(resultc <- MortalityLaw(x=x,mx=m,law="vandermaen"))
p1 <- as.numeric(coef(resultc)[1])-(min(x)-1)*as.numeric(coef(resultc)[2])+(min(x)-1)^2*as.numeric(coef(resultc)[3])
p2 <- as.numeric(coef(resultc)[2])-2*(min(x)-1)*as.numeric(coef(resultc)[3])
p3 <- as.numeric(coef(resultc)[3])
p4 <- as.numeric(coef(resultc)[4])
p5 <- as.numeric(coef(resultc)[5])+min(x)-1
error <- f(c(p1,p2,p3,p4,p5))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
oc = ifelse (is.finite(error)&&p1>=0&&p5>=0,error,Inf)
if (all(!is.finite(c(oa,ob,oc)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
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
A <- p1
B <- p2
C <- p3
I <- p4
N <- p5
}
fitted <- A+B*x+C*x^2+I/(N-x)
structure(
list(curve="Van der Maen",x=x,m=m,w=w,A=A,B=B,C=C,I=I,N=N,fitted=fitted,diagnostics=diagnostics),
class="cFit8"
)
}

#' @export
coef.cFit8 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,I=object$I,N=object$N)
}

#' @export
fitted.cFit8 <- .fittedFit

#' @export
predict.cFit8 <- function(object,newdata,...) {
object$A+object$B*newdata+object$C*newdata^2+object$I/(object$N-newdata)
}

#' @export
plot.cFit8 <- .plotFit

#' @export
deviance.cFit8 <- .devianceFit

#' @export
residuals.cFit8 <- .residualsFit
