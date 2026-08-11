#' @importFrom minpack.lm nls.lm
#' @importFrom stats nlminb
#' @importFrom graphics lines legend

cfit5 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(1/p[1]^((p[2]*x)^p[4])/p[2]+1/p[1]^((p[3]-x)^p[4])))^2) }
suppressWarnings(resulta <- nlminb(c(1,100,100,1),f,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
h <- function(p) { sqrt(w)*(log(m)-log(1/p[1]^((p[2]*x)^p[4])/p[2]+1/p[1]^((p[3]-x)^p[4]))) }
suppressWarnings(resultb <- nls.lm(c(1,100,100,1),h,lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf)))
oa = ifelse (is.finite(resulta$objective),resulta$objective,Inf)
ob = ifelse (is.finite(sum(resultb$fvec^2)),sum(resultb$fvec^2),Inf)
if (all(!is.finite(c(oa,ob)))) stop("all optimisation attempts are unsuccessful")
diagnostics = list(port=resulta,levenberg=resultb)
ind = which.min(c(oa,ob))
if (ind==1) {
A <- resulta$par[1]
B <- resulta$par[2]
M <- resulta$par[3]
N <- resulta$par[4]
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
M <- resultb$par[3]
N <- resultb$par[4]
}
fitted <- 1/A^((B*x)^N)/B+1/A^((M-x)^N)
structure(
list(curve="Wittstein Bumsted",x=x,m=m,w=w,A=A,B=B,M=M,N=N,fitted=fitted,diagnostics=diagnostics),
class="cFit5"
)
}

#' @export
coef.cFit5 <- function(object,...) {
c(A=object$A,B=object$B,M=object$M,N=object$N)
}

#' @export
fitted.cFit5 <- .fittedFit

#' @export
predict.cFit5 <- function(object,newdata,...) {
1/object$A^((object$B*newdata)^object$N)/object$B+1/object$A^((object$M-newdata)^object$N)
}

#' @export
plot.cFit5 <- .plotFit

#' @export
deviance.cFit5 <- .devianceFit

#' @export
residuals.cFit5 <- .residualsFit
