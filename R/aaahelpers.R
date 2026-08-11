.curves <- c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher","gompertz2","makeham2","oppermann2","thiele2","wittsteinbumsted2","perks2","weibull2","vandermaen2","beard2","heligmanpollard2","rogersplanck2","siler2","martinelle2","thatcher2")

.curves2 <- c("gompertz","makeham","perks","weibull","beard","martinelle","thatcher","gompertz2","makeham2","perks2","weibull2","beard2","martinelle2","thatcher2")

.checkinput <- function(x,m,w) {
if (anyNA(x)||anyNA(m)||anyNA(w)) { stop("x, m, and w must not contain missing values") }
if (!is.numeric(x)||!is.numeric(m)) { stop("x and m must be numeric") }
if (!is.vector(x)||!is.vector(m)) { stop("x and m must be vectors") }
if (length(x)!=length(m)) { stop("x and m must have the same length") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(m<=0)) { stop("m must be positive") }
if (length(w)!=length(x)) { stop("w must have the same length as x and m") }
if (any(w<0)) { stop("w must be non-negative") }
}

.fittedFit <- function(object,...) {
object$fitted
}

.plotFit <- function(x,...) {
plot(x$x,log(x$m),xlab="age",ylab="log death rate",pch=16,cex=0.5,bty="n")
lines(x$x,log(x$fitted))
legend("bottomright",legend=c("observed","fitted"),pch=c(16,NA),lty=c(NA,1),pt.cex=0.5,cex=0.8,bty="n")
}

.devianceFit <- function(object,...) {
sum(object$w*(log(object$m)-log(object$fitted))^2)
}

.residualsFit <- function(object,...) {
log(object$m)-log(object$fitted)
}

.checkinput2 <- function(x,M,h) {
if (anyNA(x)||anyNA(M)) { stop("x and M must not contain missing values") }
if (!is.numeric(x)||!is.numeric(M)) { stop("x and M must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M)) stop("M must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M)) stop("the number of ages must match the number of columns of M")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(M<=0)) { stop("all M values must be positive") }
if (nrow(M)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)) { stop("h must be numeric") }
if (h<1) { stop("h must be at least 1") }
}

.forecastModel <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { object$smoothforecast }
else if (which==2) { object$forecast }
}

.residualsModel <- function(object,...) {
object$standardresiduals
}

.checkinput3 <- function(x,M1,M2,h,jumpoff) {
if (anyNA(x)||anyNA(M1)||anyNA(M2)) { stop("x, M1, and M2 must not contain missing values") }
if (!is.numeric(x)||!is.numeric(M1)||!is.numeric(M2)) { stop("x and M1 and M2 must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(M1)||!is.matrix(M2)) stop("M1 and M2 must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(M1)||length(x)!=ncol(M2)) stop("the number of ages must match the number of columns of M1 and M2")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(M1<=0)||any(M2<=0)) { stop("all M1 and M2 values must be positive") }
if (nrow(M1)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
}

.forecastCModel <- function(object,which=1,...) {
if (length(which)!=1||!(which%in%c(1,2))) { stop("which must be 1 or 2") }
if (which==1) { list(smoothforecast1=object$smoothforecast1,smoothforecast2=object$smoothforecast2) }
else if (which==2) { list(forecast1=object$forecast1,forecast2=object$forecast2) }
}

.residualsCModel <- function(object,...) {
list(standardresiduals1=object$standardresiduals1,standardresiduals2=object$standardresiduals2)
}

.checkPN <- function(x,D,E,h,jumpoff) {
if (anyNA(x)||anyNA(D)||anyNA(E)) { stop("x, D, and E must not contain missing values") }
if (!is.numeric(x)||!is.numeric(D)||!is.numeric(E)) { stop("x and D and E must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(D)||!is.matrix(E)||nrow(D)!=nrow(E)) stop("D and E must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(D)||length(x)!=ncol(E)) stop("the number of ages must match the number of columns of D and E")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(D<=0)||any(E<=0)) { stop("all D and E values must be positive") }
if (nrow(D)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
}

.checkCPN <- function(x,D1,D2,E1,E2,h,jumpoff) {
if (anyNA(x)||anyNA(D1)||anyNA(D2)||anyNA(E1)||anyNA(E2)) { stop("x, D1, D2, E1, and E2 must not contain missing values") }
if (!is.numeric(x)||!is.numeric(D1)||!is.numeric(D2)||!is.numeric(E1)||!is.numeric(E2)) { stop("x, D1, D2, E1, and E2 must be numeric") }
if (!is.vector(x)) { stop("x must be a vector") }
if (!is.matrix(D1)||!is.matrix(D2)||!is.matrix(E1)||!is.matrix(E2)) stop("D1, D2, E1, and E2 must be a matrix with its rows as years and columns as ages")
if (length(x)!=ncol(D1)||length(x)!=ncol(D2)||length(x)!=ncol(E1)||length(x)!=ncol(E2)) stop("the number of ages must match the number of columns of D1, D2, E1, and E2")
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
if (any(D1<=0)||any(D2<=0)||any(E1<=0)||any(E2<=0)) { stop("all D1, D2, E1, and E2 values must be positive") }
if (nrow(D1)<20) stop("it requires at least 20 years of data for this forecast")
if (!is.numeric(h)||!is.numeric(jumpoff)) { stop("h and jumpoff must be numeric") }
if (h<1) { stop("h must be at least 1") }
if (jumpoff!=1&&jumpoff!=2) { stop("jump-off must be either 1 or 2") }
}
