#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit6 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log((p[1]+p[2]*p[3]^x)/(1+p[4]*p[3]^x)))^2) }
suppressWarnings(resulta <- nlminb(c(0.001,0.00001,1.1,0.00001),f))
suppressWarnings(resultb <- optim(c(0.001,0.00001,1.1,0.00001),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log((p[1]+p[2]*p[3]^x)/(1+p[4]*p[3]^x))) }
suppressWarnings(resultc <- nls.lm(c(0.001,0.00001,1.1,0.00001),h,lower=c(-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf)))
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
D <- resulta$par[4]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
D <- resultb$par[4]
} else if (ind==3) {
A <- resultc$par[1]
B <- resultc$par[2]
C <- resultc$par[3]
D <- resultc$par[4]
}
fitted <- (A+B*C^x)/(1+D*C^x)
structure(
list(curve="Perks",x=x,m=m,w=w,A=A,B=B,C=C,D=D,fitted=fitted,diagnostics=diagnostics),
class="Fit6"
)
}

#' @export
coef.Fit6 <- function(object,...) {
c(A=object$A,B=object$B,C=object$C,D=object$D)
}

#' @export
fitted.Fit6 <- .fittedFit

#' @export
predict.Fit6 <- function(object,newdata,...) {
(object$A+object$B*object$C^newdata)/(1+object$D*object$C^newdata)
}

#' @export
plot.Fit6 <- .plotFit

#' @export
deviance.Fit6 <- .devianceFit

#' @export
residuals.Fit6 <- .residualsFit
