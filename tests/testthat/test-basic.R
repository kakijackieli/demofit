library(testthat)
library(demofit)

test_that("MC runs",{
x <- 60:89
set.seed(123); m <- 0.0000082*exp(0.10771*c(60:89)+rnorm(30,0,0.1))
fit <- MC(x=x,m=m,curve="gompertz")
expect_true(!is.null(fit))
})

test_that("FCS runs",{
x <- 60:89
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962,
-3.8885,-3.7896,-3.6853,-3.5737,-3.4728,-3.3718,-3.2586,-3.1474,-3.0371,-2.9206,
-2.7998,-2.6845,-2.5653,-2.4581,-2.3367,-2.2159,-2.1017,-1.9941,-1.8821, -1.7697)
b <- c(0.0283,0.0321,0.0335,0.0336,0.0341,0.0358,0.0368,0.0403,0.0392,0.0395,
0.0396,0.0399,0.0397,0.0386,0.039,0.0375,0.0367,0.0368,0.035,0.0354,
0.0336,0.0323,0.0313,0.0295,0.0282,0.0265,0.024,0.0226,0.0219,0.0183)
k <- c(12.11,10.69,11.18,9.64,9.35,8.21,6.89,5.74,4.56,3.6,
3.27,2.04,1.11,-0.44,-1.05,-1.03,-1.84,-2.9,-4.03,-4.12,
-5.18,-5.64,-6,-6.51,-6.91,-6.9,-8.32,-8.53,-9.69,-9.31)
set.seed(123)
M <- exp(outer(k,b)+matrix(a,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.035))
fit <- FCS(x=x,M=M,model="LC",curve="makeham",h=30,jumpoff=2)
expect_true(!is.null(fit))
})

test_that("Gompertz fits",{
x <- 60:89
m <- 0.0000082*exp(0.10771*x)
fit <- MC(x=x,m=m,curve="gompertz")
cfit <- MC(x=x,m=m,curve="gompertz2")
expect_equal(unname(coef(fit)[1]),0.0000082,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.10771,tolerance=1e-4)
B <- coef(fit)[1]
C <- coef(fit)[2]
fv <- B*exp(C*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0000082,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.10771,tolerance=1e-4)
B <- coef(cfit)[1]
C <- coef(cfit)[2]
fv <- B*exp(C*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Makeham fits",{
x <- 60:89
m <- 0.005+0.0000082*exp(0.10771*x)
fit <- MC(x=x,m=m,curve="makeham")
cfit <- MC(x=x,m=m,curve="makeham2")
expect_equal(unname(coef(fit)[1]),0.005,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.0000082,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.10771,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
fv <- A+B*exp(C*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.005,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.0000082,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.10771,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
fv <- A+B*exp(C*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Oppermann fits",{
x <- 0:19
m <- 0.0038/sqrt(x+1)-0.0027+0.0005*sqrt(x+1)
fit <- MC(x=x,m=m,curve="oppermann")
cfit <- MC(x=x,m=m,curve="oppermann2")
expect_equal(unname(coef(fit)[1]),0.0038,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),-0.0027,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.0005,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
fv <- A/sqrt(x+1)+B+C*sqrt(x+1)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0038,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),-0.0027,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.0005,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
fv <- A/sqrt(x+1)+B+C*sqrt(x+1)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Thiele fits",{
x <- 0:89
m <- 0.0028/exp(1.5436*x)+0.0005/exp(0.5*0.0845*(x-21.7874)^2)+0.000058*exp(0.0872*x)
fit <- MC(x=x,m=m,curve="thiele")
cfit <- MC(x=x,m=m,curve="thiele2")
expect_equal(unname(coef(fit)[1]),0.0028,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),1.5436,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.0005,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.0845,tolerance=1e-4)
expect_equal(unname(coef(fit)[5]),21.7874,tolerance=1e-4)
expect_equal(unname(coef(fit)[6]),0.000058,tolerance=1e-4)
expect_equal(unname(coef(fit)[7]),0.0872,tolerance=1e-4)
A1 <- coef(fit)[1]
B1 <- coef(fit)[2]
A2 <- coef(fit)[3]
B2 <- coef(fit)[4]
C <- coef(fit)[5]
A3 <- coef(fit)[6]
B3 <- coef(fit)[7]
fv <- A1/exp(B1*x)+A2/exp(0.5*B2*(x-C)^2)+A3*exp(B3*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0028,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),1.5436,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.0005,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.0845,tolerance=1e-4)
expect_equal(unname(coef(cfit)[5]),21.7874,tolerance=1e-4)
expect_equal(unname(coef(cfit)[6]),0.000058,tolerance=1e-4)
expect_equal(unname(coef(cfit)[7]),0.0872,tolerance=1e-4)
A1 <- coef(cfit)[1]
B1 <- coef(cfit)[2]
A2 <- coef(cfit)[3]
B2 <- coef(cfit)[4]
C <- coef(cfit)[5]
A3 <- coef(cfit)[6]
B3 <- coef(cfit)[7]
fv <- A1/exp(B1*x)+A2/exp(0.5*B2*(x-C)^2)+A3*exp(B3*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Wittstein Bumsted fits",{
x <- 0:89
m <- 1/1.1853^((245.2946*x)^0.8563)/245.2946+1/1.1853^((106.2765-x)^0.8563)
fit <- MC(x=x,m=m,curve="wittsteinbumsted")
cfit <- MC(x=x,m=m,curve="wittsteinbumsted2")
expect_equal(unname(coef(fit)[1]),1.1853,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),245.2946,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),106.2765,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.8563,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
M <- coef(fit)[3]
N <- coef(fit)[4]
fv <- 1/A^((B*x)^N)/B+1/A^((M-x)^N)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),1.1853,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),245.2946,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),106.2765,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.8563,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
M <- coef(cfit)[3]
N <- coef(cfit)[4]
fv <- 1/A^((B*x)^N)/B+1/A^((M-x)^N)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Perks fits",{
x <- 20:89
m <- (0.000278+0.000062*1.0879^x)/(1+0.000163*1.0879^x)
fit <- MC(x=x,m=m,curve="perks")
cfit <- MC(x=x,m=m,curve="perks2")
expect_equal(unname(coef(fit)[1]),0.000278,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.000062,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),1.0879,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.000163,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
D <- coef(fit)[4]
fv <- (A+B*C^x)/(1+D*C^x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.000278,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.000062,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),1.0879,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.000163,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
D <- coef(cfit)[4]
fv <- (A+B*C^x)/(1+D*C^x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Weibull fits",{
x <- 50:69
m <- 0.0000000014*x^3.8843
fit <- MC(x=x,m=m,curve="weibull")
cfit <- MC(x=x,m=m,curve="weibull2")
expect_equal(unname(coef(fit)[1]),0.0000000014,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),3.8843,tolerance=1e-4)
B <- coef(fit)[1]
C <- coef(fit)[2]
fv <- B*x^C
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0000000014,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),3.8843,tolerance=1e-4)
B <- coef(cfit)[1]
C <- coef(cfit)[2]
fv <- B*x^C
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Van der Maen fits",{
x <- 80:99
m <- 3.096-0.0812*x+0.000511*x^2+7.545/(118.616-x)
fit <- MC(x=x,m=m,curve="vandermaen")
cfit <- MC(x=x,m=m,curve="vandermaen2")
expect_equal(unname(coef(fit)[1]),3.096,tolerance=1e-3)
expect_equal(unname(coef(fit)[2]),-0.0812,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.000511,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),7.545,tolerance=1e-3)
expect_equal(unname(coef(fit)[5]),118.616,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
I <- coef(fit)[4]
N <- coef(fit)[5]
fv <- A+B*x+C*x^2+I/(N-x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),3.096,tolerance=1e-3)
expect_equal(unname(coef(cfit)[2]),-0.0812,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.000511,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),7.545,tolerance=1e-3)
expect_equal(unname(coef(cfit)[5]),118.616,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
I <- coef(cfit)[4]
N <- coef(cfit)[5]
fv <- A+B*x+C*x^2+I/(N-x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Beard fits",{
x <- 70:99
m <- 0.000012*exp(0.1071*x)/(1+0.000002*exp(0.1071*x))
fit <- MC(x=x,m=m,curve="beard")
cfit <- MC(x=x,m=m,curve="beard2")
expect_equal(unname(coef(fit)[1]),0.000012,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.000002,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.1071,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
fv <- A*exp(C*x)/(1+B*exp(C*x))
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.000012,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.000002,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.1071,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
fv <- A*exp(C*x)/(1+B*exp(C*x))
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Heligman Pollard fits",{
x <- 0:89
xx <- ifelse(x==0,1e-10,x)
m <- 0.0004^((x+0.12)^0.16)+0.0006/exp(14.2*(log(x)-log(21.7))^2)+0.000049*1.1^x/(1+0.000049*1.1^x)
fit <- MC(x=x,m=m,curve="heligmanpollard")
cfit <- MC(x=x,m=m,curve="heligmanpollard2")
expect_equal(unname(coef(fit)[1]),0.0004,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.12,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.16,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.0006,tolerance=1e-4)
expect_equal(unname(coef(fit)[5]),14.2,tolerance=1e-4)
expect_equal(unname(coef(fit)[6]),21.7,tolerance=1e-4)
expect_equal(unname(coef(fit)[7]),0.000049,tolerance=1e-4)
expect_equal(unname(coef(fit)[8]),1.1,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
D <- coef(fit)[4]
E <- coef(fit)[5]
F <- coef(fit)[6]
G <- coef(fit)[7]
H <- coef(fit)[8]
fv <- A^((x+B)^C)+D/exp(E*(log(xx)-log(F))^2)+G*H^x/(1+G*H^x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0004,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.12,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.16,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.0006,tolerance=1e-4)
expect_equal(unname(coef(cfit)[5]),14.2,tolerance=1e-4)
expect_equal(unname(coef(cfit)[6]),21.7,tolerance=1e-4)
expect_equal(unname(coef(cfit)[7]),0.000049,tolerance=1e-4)
expect_equal(unname(coef(cfit)[8]),1.1,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
D <- coef(cfit)[4]
E <- coef(cfit)[5]
F <- coef(cfit)[6]
G <- coef(cfit)[7]
H <- coef(cfit)[8]
fv <- A^((x+B)^C)+D/exp(E*(log(xx)-log(F))^2)+G*H^x/(1+G*H^x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Rogers Planck fits",{
x <- 0:89
m <- 0.000147+0.0083/exp(2.27*x)+0.0027/exp(0.22*(x-20)+1/exp(0.28*(x-20)))+0.000077*exp(0.09*x)
fit <- MC(x=x,m=m,curve="rogersplanck")
cfit <- MC(x=x,m=m,curve="rogersplanck2")
expect_equal(unname(coef(fit)[1]),0.000147,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.0083,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.0027,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.000077,tolerance=1e-4)
expect_equal(unname(coef(fit)[5]),2.27,tolerance=1e-4)
expect_equal(unname(coef(fit)[6]),0.22,tolerance=1e-4)
expect_equal(unname(coef(fit)[7]),0.28,tolerance=1e-4)
expect_equal(unname(coef(fit)[8]),0.09,tolerance=1e-4)
expect_equal(unname(coef(fit)[9]),20,tolerance=1e-4)
A0 <- coef(fit)[1]
A1 <- coef(fit)[2]
A2 <- coef(fit)[3]
A3 <- coef(fit)[4]
A <- coef(fit)[5]
B <- coef(fit)[6]
C <- coef(fit)[7]
D <- coef(fit)[8]
U <- coef(fit)[9]
fv <- A0+A1/exp(A*x)+A2/exp(B*(x-U)+1/exp(C*(x-U)))+A3*exp(D*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.000147,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.0083,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.0027,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.000077,tolerance=1e-4)
expect_equal(unname(coef(cfit)[5]),2.27,tolerance=1e-4)
expect_equal(unname(coef(cfit)[6]),0.22,tolerance=1e-4)
expect_equal(unname(coef(cfit)[7]),0.28,tolerance=1e-4)
expect_equal(unname(coef(cfit)[8]),0.09,tolerance=1e-4)
expect_equal(unname(coef(cfit)[9]),20,tolerance=1e-4)
A0 <- coef(cfit)[1]
A1 <- coef(cfit)[2]
A2 <- coef(cfit)[3]
A3 <- coef(cfit)[4]
A <- coef(cfit)[5]
B <- coef(cfit)[6]
C <- coef(cfit)[7]
D <- coef(cfit)[8]
U <- coef(cfit)[9]
fv <- A0+A1/exp(A*x)+A2/exp(B*(x-U)+1/exp(C*(x-U)))+A3*exp(D*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Siler fits",{
x <- 0:89
m <- 0.0032/exp(2.03*x)+0.000008+0.000074*exp(0.0853*x)
fit <- MC(x=x,m=m,curve="siler")
cfit <- MC(x=x,m=m,curve="siler2")
expect_equal(unname(coef(fit)[1]),0.0032,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),2.03,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.000008,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.000074,tolerance=1e-4)
expect_equal(unname(coef(fit)[5]),0.0853,tolerance=1e-4)
A1 <- coef(fit)[1]
B1 <- coef(fit)[2]
A2 <- coef(fit)[3]
A3 <- coef(fit)[4]
B3 <- coef(fit)[5]
fv <- A1/exp(B1*x)+A2+A3*exp(B3*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0032,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),2.03,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.000008,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.000074,tolerance=1e-4)
expect_equal(unname(coef(cfit)[5]),0.0853,tolerance=1e-4)
A1 <- coef(cfit)[1]
B1 <- coef(cfit)[2]
A2 <- coef(cfit)[3]
A3 <- coef(cfit)[4]
B3 <- coef(cfit)[5]
fv <- A1/exp(B1*x)+A2+A3*exp(B3*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Martinelle fits",{
x <- 20:89
m <- (0.0048+0.000061*exp(0.0879*x))/(1+0.5596*exp(0.0879*x))+0.000057*exp(0.0879*x)
fit <- MC(x=x,m=m,curve="martinelle")
cfit <- MC(x=x,m=m,curve="martinelle2")
expect_equal(unname(coef(fit)[1]),0.0048,tolerance=1e-3)
expect_equal(unname(coef(fit)[2]),0.000061,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.0879,tolerance=1e-3)
expect_equal(unname(coef(fit)[4]),0.5596,tolerance=1e-3)
expect_equal(unname(coef(fit)[5]),0.000057,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
D <- coef(fit)[4]
E <- coef(fit)[5]
fv <- (A+B*exp(C*x))/(1+D*exp(C*x))+E*exp(C*x)
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.0048,tolerance=1e-3)
expect_equal(unname(coef(cfit)[2]),0.000061,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.0879,tolerance=1e-3)
expect_equal(unname(coef(cfit)[4]),0.5596,tolerance=1e-3)
expect_equal(unname(coef(cfit)[5]),0.000057,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
D <- coef(cfit)[4]
E <- coef(cfit)[5]
fv <- (A+B*exp(C*x))/(1+D*exp(C*x))+E*exp(C*x)
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("Thatcher fits",{
x <- 80:104
m <- 0.02+0.0001*exp(0.14*x)/(1+0.000001*exp(0.14*x))
fit <- MC(x=x,m=m,curve="thatcher")
cfit <- MC(x=x,m=m,curve="thatcher2")
expect_equal(unname(coef(fit)[1]),0.02,tolerance=1e-4)
expect_equal(unname(coef(fit)[2]),0.0001,tolerance=1e-4)
expect_equal(unname(coef(fit)[3]),0.14,tolerance=1e-4)
expect_equal(unname(coef(fit)[4]),0.000001,tolerance=1e-4)
A <- coef(fit)[1]
B <- coef(fit)[2]
C <- coef(fit)[3]
D <- coef(fit)[4]
fv <- A+B*exp(C*x)/(1+D*exp(C*x))
expect_equal(fitted(fit),fv,tolerance=1e-4)
expect_equal(unname(coef(cfit)[1]),0.02,tolerance=1e-4)
expect_equal(unname(coef(cfit)[2]),0.0001,tolerance=1e-4)
expect_equal(unname(coef(cfit)[3]),0.14,tolerance=1e-4)
expect_equal(unname(coef(cfit)[4]),0.000001,tolerance=1e-4)
A <- coef(cfit)[1]
B <- coef(cfit)[2]
C <- coef(cfit)[3]
D <- coef(cfit)[4]
fv <- A+B*exp(C*x)/(1+D*exp(C*x))
expect_equal(fitted(cfit),fv,tolerance=1e-4)
})

test_that("LCS fits",{
x <- 60:69
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
b <- c(0.0801,0.0909,0.0948,0.0951,0.0965,0.1014,0.1042,0.1141,0.1111,0.1118)
k <- c(8.46,7.04,7.53,5.99,5.70,4.56,3.24,2.09,0.91,-0.05,
-0.38,-1.61,-2.54,-4.09,-4.70,-4.68,-5.49,-6.55,-7.66,-7.77)
M <- exp(outer(k,b)+matrix(a,nrow=20,ncol=10,byrow=TRUE))
fit <- LCS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$alpha,a,tolerance=1e-4)
expect_equal(coef(fit)$beta,b,tolerance=1e-4)
expect_equal(coef(fit)$kappa,k,tolerance=1e-4)
})

test_that("RHS fits",{
x <- 60:69
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
b <- c(0.0801,0.0909,0.0948,0.0951,0.0965,0.1014,0.1042,0.1141,0.1111,0.1118)
k <- c(8.46,7.04,7.53,5.99,5.70,4.56,3.24,2.09,0.91,-0.05,
-0.38,-1.61,-2.54,-4.09,-4.70,-4.68,-5.49,-6.55,-7.66,-7.77)
g <- c(0,0,0,0,0,0.05,-0.05,0,0,0,0,0,0,0,0,0.06,-0.06,0,0,0,0,0,0,0,0,0,0,0,0)
M <- outer(k,b)+matrix(a,nrow=20,ncol=10,byrow=TRUE)
M <- M+g[row(M)-col(M)+10]
M <- exp(M)
fit <- RHS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$alpha,a,tolerance=1e-4)
expect_equal(coef(fit)$beta,b,tolerance=1e-4)
expect_equal(coef(fit)$kappa,k,tolerance=1e-3)
expect_equal(coef(fit)$gamma,g,tolerance=1e-2)
})

test_that("APCS fits",{
x <- 60:69
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
k <- c(8.46,7.04,7.53,5.99,5.70,4.56,3.24,2.09,0.91,-0.05,
-0.38,-1.61,-2.54,-4.09,-4.70,-4.68,-5.49,-6.55,-7.66,-7.77)
g <- c(0,0,0,0,0,0.05,-0.05,0,0,0,0,0,0,0,0,0.06,-0.06,0,0,0,0,0,0,0,0,0,0,0,0)
M <- matrix(a,nrow=20,ncol=10,byrow=TRUE)+matrix(k,nrow=20,ncol=10,byrow=FALSE)
M <- M+g[row(M)-col(M)+10]
M <- exp(M)
fit <- APCS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$alpha,a,tolerance=1e-4)
expect_equal(coef(fit)$kappa,k,tolerance=1e-4)
expect_equal(coef(fit)$gamma,g,tolerance=1e-3)
})

test_that("CBDS fits",{
x <- 60:69
k1 <- -2.97-0.0245*(0:19)
k2 <- 0.101+0.000345*(0:19)
M <- exp(matrix(k1,nrow=20,ncol=10,byrow=FALSE)+outer(k2,(x-mean(x))))
fit <- CBDS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$kappa1,k1,tolerance=1e-4)
expect_equal(coef(fit)$kappa2,k2,tolerance=1e-4)
})

test_that("CBDCS fits",{
x <- 60:69
k1 <- -2.97-0.0245*(0:19)
k2 <- 0.101+0.000345*(0:19)
g <- c(0,0,0,0,0,0.05,-0.05,0,0,0,0,0,0,0,0,0.06,-0.06,0,0,0,0,0,0,0,0,0,0,0,0)
M <- matrix(k1,nrow=20,ncol=10,byrow=FALSE)+outer(k2,(x-mean(x)))
M <- M+g[row(M)-col(M)+10]
M <- exp(M)
fit <- CBDCS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$kappa1,k1,tolerance=1e-4)
expect_equal(coef(fit)$kappa2,k2,tolerance=1e-4)
expect_equal(coef(fit)$gamma,g,tolerance=1e-3)
})

test_that("CBDQCS fits",{
x <- 60:69
s2 <- mean((x-mean(x))^2)
k1 <- -2.97-0.0245*(0:19)
k2 <- 0.101+0.000345*(0:19)
k3 <- (0:19)/1000
g <- c(0,0,0,0,0,0.05,-0.05,0,0,0,0,0,0,0,0,0.06,-0.06,0,0,0,0,0,0,0,0,0,0,0,0)
M <- matrix(k1,nrow=20,ncol=10,byrow=FALSE)+outer(k2,(x-mean(x)))+outer(k3,(x-mean(x))^2-s2)
M <- M+g[row(M)-col(M)+10]
M <- exp(M)
fit <- CBDQCS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$kappa1,k1,tolerance=1e-3)
expect_equal(coef(fit)$kappa2,k2,tolerance=1e-3)
expect_equal(coef(fit)$kappa3,k3,tolerance=1e-4)
expect_equal(coef(fit)$gamma,g,tolerance=1e-2)
})

test_that("STARS fits",{
x <- 60:69
mu <- c(-0.0226,0.0254,0.0178,0.0616,0.0606,0.0585,0.0183,0.0601,0.1040,0.0538)
R <- matrix(0,nrow=10,ncol=10)
R[1,1] <- 1
R[2,1] <- 0.706; R[2,2] <- 0.294
R[3,1] <- 0.059; R[3,2] <- 0.424; R[3,3] <- 0.517
R[4,2] <- 0.335; R[4,3] <- 0.415; R[4,4] <- 0.250
R[5,3] <- 0.122; R[5,4] <- 0.830; R[5,5] <- 0.048
R[6,4] <- 0.125; R[6,5] <- 0.860; R[6,6] <- 0.015
R[7,5] <- 0.231; R[7,6] <- 0.093; R[7,7] <- 0.676
R[8,6] <- 0.065; R[8,7] <- 0.910; R[8,8] <- 0.025
R[9,7] <- 0.591; R[9,8] <- 0.384; R[9,9] <- 0.025
R[10,8] <- 0.010; R[10,9] <- 0.980; R[10,10] <- 0.010
M <- matrix(0,nrow=20,ncol=10)
M[1,] <- c(0.0165,0.0174,0.0190,0.0202,0.0219,0.0236,0.0250,0.0279,0.0289,0.0315)
for (t in 2:20) { M[t,] <- exp(mu+R%*%log(M[t-1,])) }
fit <- STARS(x=x,M=M,curve="makeham")
expect_equal(coef(fit)$mu,mu,tolerance=1e-4)
expect_equal(coef(fit)$R,R,tolerance=1e-4)
})

test_that("edge cases",{
x <- 60:89
m <- rep(0,30)
M <- matrix(0,nrow=30,ncol=30)
expect_error(MC(x=x,m=m,curve="gompertz"))
expect_error(FCS(x=x,M=M,curve="gompertz"))
M <- matrix(0.01,nrow=10,ncol=30)
expect_error(FCS(x=x,M=M,curve="gompertz"))
m <- rep(NA,30)
M <- matrix(NA,nrow=30,ncol=30)
expect_error(MC(x=x,m=m,curve="gompertz"))
expect_error(FCS(x=x,M=M,curve="gompertz"))
})

test_that("valid forecasts",{
x <- 60:69
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
b <- c(0.0801,0.0909,0.0948,0.0951,0.0965,0.1014,0.1042,0.1141,0.1111,0.1118)
k <- c(8.46,7.04,7.53,5.99,5.70,4.56,3.24,2.09,0.91,-0.05,
-0.38,-1.61,-2.54,-4.09,-4.70,-4.68,-5.49,-6.55,-7.66,-7.77)
M <- exp(outer(k,b)+matrix(a,nrow=20,ncol=10,byrow=TRUE))
fit <- FCS(x=x,M=M,model="LC",curve="makeham")
expect_false(any(is.na(forecast::forecast(fit))))
expect_true(all(is.finite(forecast::forecast(fit))))
expect_true(all(forecast::forecast(fit)>0))
})

test_that("ENS runs",{
x <- 60:69
a <- c(-4.8499,-4.7676,-4.6719,-4.5722,-4.4847,-4.3841,-4.2813,-4.1863,-4.0861,-3.9962)
b <- c(0.0801,0.0909,0.0948,0.0951,0.0965,0.1014,0.1042,0.1141,0.1110,0.1118)
k <- c(12.11,10.69,11.18,9.64,9.35,8.21,6.89,5.74,4.56,3.60,
3.27,2.04,1.11,-0.44,-1.05,-1.03,-1.84,-2.90,-4.03,-4.12,
-5.18,-5.64,-6.00,-6.51,-6.91,-6.90,-8.32,-8.53,-9.69,-9.31)
set.seed(123)
M <- exp(outer(k,b)+matrix(a,nrow=30,ncol=10,byrow=TRUE)+rnorm(300,0,0.035))
fit <- ENS(x=x,M=M,curve="makeham",h=30,jumpoff=2)
expect_true(!is.null(fit))
expect_false(any(is.na(forecast::forecast(fit))))
expect_true(all(is.finite(forecast::forecast(fit))))
expect_true(all(forecast::forecast(fit)>0))
})

test_that("CFMS runs",{
x <- 60:89
a1 <- c(-5.18,-5.12,-4.98,-4.92,-4.82,-4.73,-4.66,-4.53,-4.45,-4.35,
-4.26,-4.17,-4.05,-3.95,-3.84,-3.73,-3.65,-3.52,-3.40,-3.29,
-3.14,-3.02,-2.88,-2.76,-2.64,-2.49,-2.37,-2.25,-2.12,-2.00)
a2 <- c(-4.78,-4.68,-4.57,-4.49,-4.39,-4.29,-4.19,-4.10,-4.00,-3.89,
-3.80,-3.69,-3.60,-3.49,-3.39,-3.29,-3.17,-3.07,-2.96,-2.85,
-2.71,-2.62,-2.49,-2.37,-2.26,-2.14,-2.04,-1.91,-1.82,-1.72)
B <- c(0.0381,0.0340,0.0420,0.0389,0.0423,0.0414,0.0406,0.0393,0.0415,0.0400,
0.0411,0.0362,0.0387,0.0381,0.0384,0.0385,0.0356,0.0314,0.0317,0.0337,
0.0316,0.0298,0.0284,0.0270,0.0248,0.0262,0.0205,0.0215,0.0142,0.0145)
K <- c(9.66,9.89,10.66,9.83,9.52,7.39,7.64,6.36,2.32,4.18,
2.91,-0.61,0.28,-0.38,-1.79,-3.34,-1.74,-3.50,-4.28,-4.77,
-4.98,-7.13,-5.09,-6.41,-5.56,-5.65,-6.12,-5.64,-7.35,-6.28)
b1 <- c(0.0012,-0.0033,0.0523,0.0161,0.0529,0.0220,0.0312,0.0437,0.0709,0.0444,
0.0398,0.0361,0.0403,0.0396,0.0506,0.0315,0.0428,0.0261,0.0384,0.0388,
0.0300,0.0269,0.0275,0.0256,0.0239,0.0421,0.0314,0.0284,0.0174,0.0314)
k1 <- c(-1.24,-1.38,-3.48,-2.51,-1.32,-1.90,-3.42,-0.94,0.24,-0.48,
-0.26,2.70,1.39,-0.46,1.74,2.53,0.90,1.43,0.76,2.48,
0.74,2.32,0.42,1.69,-0.64,1.30,0.19,-0.69,-1.11,-1.01)
b2 <- c(-0.0014,0.0272,0.0083,0.0273,0.0209,0.0253,0.0144,0.0333,0.0460,0.0439,
0.0439,0.0674,0.0331,0.0443,0.0312,0.0240,0.0570,0.0312,0.0403,0.0376,
0.0500,0.0289,0.0466,0.0418,0.0349,0.0149,0.0366,0.0178,0.0361,0.0372)
k2 <- c(2.35,0.62,-0.38,0.12,0.00,0.80,-1.39,0.38,2.47,0.40,
0.76,3.06,1.42,-0.73,0.79,1.94,0.12,0.60,-0.43,0.29,
0.17,0.98,-1.01,-0.13,-2.46,-1.24,-1.65,-2.48,-2.32,-3.06)
set.seed(123)
M1 <- exp(outer(k1,b1)+outer(K,B)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
M2 <- exp(outer(k2,b2)+outer(K,B)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
fit <- CFMS(x=x,M1=M1,M2=M2,curve="makeham",h=30,jumpoff=2)
expect_true(!is.null(fit))
expect_false(any(is.na(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(is.finite(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(forecast::forecast(fit)$smoothforecast2>0))
})

test_that("CFM2S runs",{
x <- 60:89
a1 <- c(-5.18,-5.12,-4.98,-4.92,-4.82,-4.73,-4.66,-4.53,-4.45,-4.35,
-4.26,-4.17,-4.05,-3.95,-3.84,-3.73,-3.65,-3.52,-3.40,-3.29,
-3.14,-3.02,-2.88,-2.76,-2.64,-2.49,-2.37,-2.25,-2.12,-2.00)
a2 <- c(-4.78,-4.68,-4.57,-4.49,-4.39,-4.29,-4.19,-4.10,-4.00,-3.89,
-3.80,-3.69,-3.60,-3.49,-3.39,-3.29,-3.17,-3.07,-2.96,-2.85,
-2.71,-2.62,-2.49,-2.37,-2.26,-2.14,-2.04,-1.91,-1.82,-1.72)
B <- c(0.0381,0.0340,0.0420,0.0389,0.0423,0.0414,0.0406,0.0393,0.0415,0.0400,
0.0411,0.0362,0.0387,0.0381,0.0384,0.0385,0.0356,0.0314,0.0317,0.0337,
0.0316,0.0298,0.0284,0.0270,0.0248,0.0262,0.0205,0.0215,0.0142,0.0145)
K <- c(9.66,9.89,10.66,9.83,9.52,7.39,7.64,6.36,2.32,4.18,
2.91,-0.61,0.28,-0.38,-1.79,-3.34,-1.74,-3.50,-4.28,-4.77,
-4.98,-7.13,-5.09,-6.41,-5.56,-5.65,-6.12,-5.64,-7.35,-6.28)
b <- c(-0.00010,0.01195,0.03030,0.02170,0.03690,0.02365,0.02280,0.03850,0.05845,0.04415,
0.04185,0.05175,0.03670,0.04195,0.04090,0.02775,0.04990,0.02865,0.03935,0.03820,
0.04000,0.02790,0.03705,0.03370,0.02940,0.02850,0.03400,0.02310,0.02675,0.03430)
k1 <- c(-1.24,-1.38,-3.48,-2.51,-1.32,-1.90,-3.42,-0.94,0.24,-0.48,
-0.26,2.70,1.39,-0.46,1.74,2.53,0.90,1.43,0.76,2.48,
0.74,2.32,0.42,1.69,-0.64,1.30,0.19,-0.69,-1.11,-1.01)
k2 <- c(2.35,0.62,-0.38,0.12,0.00,0.80,-1.39,0.38,2.47,0.40,
0.76,3.06,1.42,-0.73,0.79,1.94,0.12,0.60,-0.43,0.29,
0.17,0.98,-1.01,-0.13,-2.46,-1.24,-1.65,-2.48,-2.32,-3.06)
set.seed(123)
M1 <- exp(outer(k1,b)+outer(K,B)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
M2 <- exp(outer(k2,b)+outer(K,B)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
fit <- CFM2S(x=x,M1=M1,M2=M2,curve="makeham",h=30,jumpoff=2)
expect_true(!is.null(fit))
expect_false(any(is.na(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(is.finite(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(forecast::forecast(fit)$smoothforecast2>0))
})

test_that("CAES runs",{
x <- 60:89
a1 <- c(-5.18,-5.12,-4.98,-4.92,-4.82,-4.73,-4.66,-4.53,-4.45,-4.35,
-4.26,-4.17,-4.05,-3.95,-3.84,-3.73,-3.65,-3.52,-3.40,-3.29,
-3.14,-3.02,-2.88,-2.76,-2.64,-2.49,-2.37,-2.25,-2.12,-2.00)
a2 <- c(-4.78,-4.68,-4.57,-4.49,-4.39,-4.29,-4.19,-4.10,-4.00,-3.89,
-3.80,-3.69,-3.60,-3.49,-3.39,-3.29,-3.17,-3.07,-2.96,-2.85,
-2.71,-2.62,-2.49,-2.37,-2.26,-2.14,-2.04,-1.91,-1.82,-1.72)
b <- c(0.0381,0.0340,0.0420,0.0389,0.0423,0.0414,0.0406,0.0393,0.0415,0.0400,
0.0411,0.0362,0.0387,0.0381,0.0384,0.0385,0.0356,0.0314,0.0317,0.0337,
0.0316,0.0298,0.0284,0.0270,0.0248,0.0262,0.0205,0.0215,0.0142,0.0145)
k1 <- c(8.68,8.34,7.99,6.87,8.18,5.73,4.83,5.20,2.74,3.22,
2.99,1.59,1.67,-0.65,-0.39,-1.07,-0.95,-2.78,-3.46,-2.45,
-4.12,-4.66,-4.98,-4.58,-6.30,-4.39,-5.56,-6.52,-8.26,-6.92)
k2 <- c(11.81,11.01,10.59,10.40,9.75,8.15,6.07,6.45,4.60,4.57,
4.15,1.49,1.77,-1.08,-1.44,-0.96,-1.66,-2.25,-4.67,-4.62,
-4.38,-6.37,-6.27,-6.91,-8.22,-7.35,-8.39,-7.87,-9.72,-8.65)
set.seed(123)
M1 <- exp(outer(k1,b)+matrix(a1,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
M2 <- exp(outer(k2,b)+matrix(a2,nrow=30,ncol=30,byrow=TRUE)+rnorm(900,0,0.07))
fit <- CAES(x=x,M1=M1,M2=M2,curve="makeham",h=30,jumpoff=2)
expect_true(!is.null(fit))
expect_false(any(is.na(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(is.finite(forecast::forecast(fit)$smoothforecast1)))
expect_true(all(forecast::forecast(fit)$smoothforecast2>0))
})
