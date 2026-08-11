#' @importFrom minpack.lm nls.lm
#' @importFrom stats optim nlminb
#' @importFrom graphics lines legend

fit3 <- function(x,m,w) {
.checkinput(x,m,w)
f <- function(p) { sum(w *(log(m)-log(p[1]/sqrt(x+1)-p[2]+p[3]*sqrt(x+1)))^2) }
suppressWarnings(resulta <- nlminb(c(0.01,0.001,0.001),f))
suppressWarnings(resultb <- optim(c(0.01,0.001,0.001),f,method="Nelder-Mead"))
h <- function(p) { sqrt(w)*(log(m)-log(p[1]/sqrt(x+1)-p[2]+p[3]*sqrt(x+1))) }
suppressWarnings(resultc <- nls.lm(c(0.01,0.001,0.001),h,lower=c(-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf)))
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
} else if (ind==2) {
A <- resultb$par[1]
B <- resultb$par[2]
C <- resultb$par[3]
} else if (ind==3) {
A <- resultc$par[1]
B <- resultc$par[2]
C <- resultc$par[3]
}
fitted <- A/sqrt(x+1)-B+C*sqrt(x+1)
structure(
list(curve="Oppermann",x=x,m=m,w=w,A=A,B=B,C=C,fitted=fitted,diagnostics=diagnostics),
class="Fit3"
)
}

#' @export
coef.Fit3 <- function(object,...) {
c(A=object$A,B=-object$B,C=object$C)
}

#' @export
fitted.Fit3 <- .fittedFit

#' @export
predict.Fit3 <- function(object,newdata,...) {
object$A/sqrt(newdata+1)-object$B+object$C*sqrt(newdata+1)
}

#' @export
plot.Fit3 <- .plotFit

#' @export
deviance.Fit3 <- .devianceFit

#' @export
residuals.Fit3 <- .residualsFit
