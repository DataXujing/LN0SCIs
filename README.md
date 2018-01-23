
[![CRAN](https://img.shields.io/cran/v/devtools.svg)]( https://CRAN.R-project.org/package=LN0SCIs )
[![PyPI](https://img.shields.io/pypi/v/nine.svg)](https://pypi.python.org/pypi/LN0SCIs)
# LN0SCIs   <a href="https://github.com/DataXujing/"><img src="pic/log.png" align="right" alt="logo" height="200" width="350" /></a>


## Introduction

This R package based on the paper of Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros by Xu et al. It provides serveral methods for construct simultaneous confidence intervals for ratios of means of log-normal populations with excess zeros. At last, we select 4 excellent methods which based on generalized pivotal quantity with order statistics and two-step MOVER intervals.
For the convenience of use, we make a R package called `LN0SCIs`, and it also has a python version package: <https://pypi.python.org/pypi/LN0SCIs>. You can use `vignettes('LN0SCIs-tutorial')` in R  to read the tototial.


## Methods

We provaide four main functions in our LN0SCIs packages, FGW(),FGH(),MOVERW() and MOVERH(), if you want to  deep understanding these four methods,you can read our paper: Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros. the code we trust in GitHub. If you want to know how to realize them,you can read the source code.


## Example

+ FGW()

```r
library(LN0SCIs)

# params setting
p<-c(0.1,0.15,0.1,0.6)
n<-c(30,15,10,50)
mu<-c(1,1.3,2,0)
sigma<-c(1,1,1,2)
C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))

N<-1000
FGW(n,p,mu,sigma,N,C2 = C2) #base function
```

```r
[1] "====================Method: FGW====================="
[1] "The Simultaneous Confidence Intervals are:          "
               [LCL,UCL]
1   [-1.235113,2.848869]
2   [-0.441577,7.030192]
3   [-3.57776,-2.108937]
4    [-1.86122,6.480295]
5  [-5.843864,-1.640372]
6 [-10.008059,-2.477454]
[1] "**********************Time**************************"
Time difference of 53.041 secs
```
+ FGH()


```r

p<-c(0.1,0.15,0.1,0.6)
n<-c(30,15,10,50)
mu<-c(1,1.3,2,0)
sigma<-c(1,1,1,2)
C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))

N<-1000;
FGH(n,p,mu,sigma,N,C2 = C2)
```

```r
[1] "====================Method: FGH====================="
[1] "The Simultaneous Confidence Intervals are:          "
             [LCL,UCL]
1  [-0.996117,1.07758]
2  [-0.36159,2.378184]
3  [-2.87464,1.273622]
4  [-0.302924,2.22018]
5 [-2.954485,1.098756]
6 [-4.086911,0.345388]
[1] "**********************Time**************************"
Time difference of 16.244 secs
```

+ MOVERW()

```r
p <- c(0.1,0.15,0.1,0.6)
n <- c(30,15,10,50)
mu <- c(1,1.3,2,0)
sigma <- c(1,1,1,2)
C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
N <- 1000

MOVERW(n,p,mu,sigma,N,C2 = C2)
```

```r
[1] "====================Method: MOVERW=================="
[1] "The Simultaneous Confidence Intervals are:          "
             [LCL,UCL]
1 [-1.351672,0.720717]
2  [-1.11909,2.371603]
3 [-1.882617,2.617392]
4  [-0.86023,2.700871]
5 [-1.618718,2.936072]
6 [-2.967078,2.541986]
[1] "**********************Time**************************"
Time difference of 42.577 secs
```

+ MOVERH()

```r
p<-c(0.1,0.15,0.1,0.6)
n<-c(30,15,10,50)
mu<-c(1,1.3,2,0)
sigma<-c(1,1,1,2)
C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))

N<-1000;
MOVERH(n,p,mu,sigma,N,C2 = C2)
```

```r
[1] "====================Method: MOVERH=================="
[1] "The Simultaneous Confidence Intervals are:          "
             [LCL,UCL]
1 [-1.334835,1.683172]
2 [-0.806956,2.953145]
3 [-2.683212,1.739062]
4 [-1.130242,2.909716]
5 [-2.979693,1.684558]
6 [-4.145004,1.117301]
[1] "**********************Time**************************"
Time difference of 11.152 secs
```



