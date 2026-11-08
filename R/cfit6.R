#' @importFrom minpack.lm nls.lm
#' @importFrom MortalityLaws MortalityLaw
#' @importFrom stats nlminb coef
#' @importFrom graphics lines legend

cfit6 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log((p[1]+p[2]*p[3]^x)/(1+p[4]*p[3]^x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,1.1,0.00001),f,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log((p[1]+p[2]*p[3]^x)/(1+p[4]*p[3]^x))) }
suppressWarnings(resultb <- nls.lm(c(0.001,0.00001,1.1,0.00001),h,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
suppressWarnings(resultc <- MortalityLaw(x=x,mx=m,law="perks"))
p1 <- as.numeric(coef(resultc)[1])
p2 <- as.numeric(coef(resultc)[2])*as.numeric(coef(resultc)[3])^(1-min(x))
p3 <- as.numeric(coef(resultc)[3])
p4 <- as.numeric(coef(resultc)[4])*as.numeric(coef(resultc)[3])^(1-min(x))
error <- f(c(p1,p2,p3,p4))
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
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
D <- resultb$par[4]
} else if (ind==3) {
A <- p1
B <- p2
C <- p3
D <- p4
}
fitted <- (A+B*C^x)/(1+D*C^x)
structure(
list(curve="Perks",x=x,m=m,w=w,A=A,B=B,C=C,D=D,fitted=fitted,diagnostics=diagnostics),
class="cFit6"
)
}

#' @export
coef.cFit6 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D)
}

#' @export
fitted.cFit6 <- .fittedFit

#' @export
predict.cFit6 <- function(object,newdata,...) {
(object$A+object$B*object$C^newdata)/(1+object$D*object$C^newdata)
}

#' @export
plot.cFit6 <- .plotFit

#' @export
deviance.cFit6 <- .devianceFit

#' @export
residuals.cFit6 <- .residualsFit
