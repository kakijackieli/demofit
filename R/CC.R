#' Mortality curve calculation
#'
#' Calculates mortality rates given mortality curve parameter values.
#'
#' @param x vector of ages.
#' @param par vector of mortality curve parameter values.
#' @param curve name of mortality curve (including gompertz, makeham, oppermann, thiele, wittsteinbumsted, perks, weibull, vandermaen, beard, heligmanpollard, rogersplanck, siler, martinelle, thatcher).
#'
#' @details
#' See \code{MC()} for more details of different mortality curves.
#'
#' @return 
#' Calculated mortality rates based on selected mortality curve and parameter values.
#'
#' @examples
#' CC(x=60:89,par=c(0.0000082,0.10771),curve="gompertz")
#'
#' @export
CC <- function(x,par,curve=c("gompertz","makeham","oppermann","thiele","wittsteinbumsted","perks","weibull","vandermaen","beard","heligmanpollard","rogersplanck","siler","martinelle","thatcher")) {
if (!is.numeric(x)||!is.numeric(par)) { stop("both x and par must be numeric") }
if (!is.vector(x)||!is.vector(par)) { stop("x and par must be vectors") }
if (is.unsorted(x,strictly=TRUE)) { stop("x must be in ascending order") }
if (any(x<0)) { stop("x must be non-negative") }
curve <- tryCatch(match.arg(curve),error = function(e) { stop("invalid curve choice") })
tryCatch({
parcount <- list(gompertz=2,makeham=3,oppermann=3,thiele=7,wittsteinbumsted=4,perks=4,weibull=2,vandermaen=5,beard=3,heligmanpollard=8,rogersplanck=9,siler=5,martinelle=5,thatcher=4)
if (length(par)!=parcount[[curve]]) { stop("the number of parameters is incorrect for the curve") }
xx <- ifelse(x==0,1e-10,x)
switch(curve,
gompertz = { B <- par[1]; C <- par[2]; B*exp(C*x) },
makeham  = { A <- par[1]; B <- par[2]; C <- par[3]; A+B*exp(C*x) },
oppermann = { A <- par[1]; B <- par[2]; C <- par[3]; A/sqrt(x+1)+B+C*sqrt(x+1) },
thiele = { A1 <- par[1]; B1 <- par[2]; A2 <- par[3]; B2 <- par[4]; C <- par[5]; A3 <- par[6]; B3 <- par[7]; A1/exp(B1*x)+A2/exp(0.5*B2*(x-C)^2)+A3*exp(B3*x) },
wittsteinbumsted = { A <- par[1]; B <- par[2]; M <- par[3]; N <- par[4]; 1/A^((B*x)^N)/B+1/A^((M-x)^N) },
perks = { A <- par[1]; B <- par[2]; C <- par[3]; D <- par[4]; (A+B*C^x)/(1+D*C^x) },
weibull = { B <- par[1]; C <- par[2]; B*x^C },
vandermaen = { A <- par[1]; B <- par[2]; C <- par[3]; I <- par[4]; N <- par[5]; A+B*x+C*x^2+I/(N-x) },
beard = { A <- par[1]; B <- par[2]; C <- par[3]; A*exp(C*x)/(1+B*exp(C*x)) },
heligmanpollard = { A <- par[1]; B <- par[2]; C <- par[3]; D <- par[4]; E <- par[5]; F <- par[6]; G <- par[7]; H <- par[8]; A^((x+B)^C)+D/exp(E*(log(xx)-log(F))^2)+G*H^x/(1+G*H^x) },
rogersplanck = { A0 <- par[1]; A1 <- par[2]; A2 <- par[3]; A3 <- par[4]; A <- par[5]; B <- par[6]; C <- par[7]; D <- par[8]; U <- par[9]; A0+A1/exp(A*x)+A2/exp(B*(x-U)+1/exp(C*(x-U)))+A3*exp(D*x) },
siler = { A1 <- par[1]; B1 <- par[2]; A2 <- par[3]; A3 <- par[4]; B3 <- par[5]; A1/exp(B1*x)+A2+A3*exp(B3*x) },
martinelle = { A <- par[1]; B <- par[2]; C <- par[3]; D <- par[4]; E <- par[5]; (A+B*exp(C*x))/(1+D*exp(C*x))+E*exp(C*x) },
thatcher = { A <- par[1]; B <- par[2]; C <- par[3]; D <- par[4]; A+B*exp(C*x)/(1+D*exp(C*x)) }
)
}, error = function(e) { stop("curve calculation is unsuccessful - please make sure the parameters and age range are suitable for the curve") })
}
