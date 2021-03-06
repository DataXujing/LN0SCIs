#' LN0SCIs
#'
#' Construct the simultaneous confidence intervals for
#' ratios of means of Log-normal populations with excess zeros. It also has a Python module that do the same thing,
#' it can be applied to multiple comparisons of parameters of any k mixture distribution. And it provide four methods,
#' the method based on generalized pivotal quantity with order statistics (\code{\link{FGH}} and \code{\link{FGW}}),
#' and the method based on two-step MOVER intervals (\code{\link{MOVERW}} and \code{\link{MOVERH}}).
#'
#'
#' @name LN0SCIs
#' @aliases LN0SCIs
#' @docType package
#' @author Jing Xu, Xinmin Li, Hua Liang
#'
#' @details At present, these four function perform better than other methods
#' that can be used to calculate the simultaneous confidence interval of log-normal populations with excess zeros.
#'
#' @import stats
#'
#' @seealso
#'
#' [1] Besag I, Green P, Higdon D, Mengersen K, 1995. Bayesian computation and Stochastic-systems.
#'
#' [2] Hannig J, Abdel-Karim A, Iyer H, 2006. Simultaneous fiducial generalized confidence
#' intervals for ratios of means of lognormal distribution.
#'
#' [3] Hannig J, Lee T C M, 2009.Generalized fiducial inference for wavelet regression.
#'
#' [4] Li X, Zhou X, Tian L, 2013. Interval estimation for the mean of lognormal data with excess zeros.
#'
#' [5] Schaarschmidt F, 2013. Simultaneous confidence intervals for multiple comparisons among expected values of log-normal variables.
#'
#' [6] Jing Xu, Xinmin Li, Hua Liang. Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros.
#'
#'@examples
#'
#'\dontrun{
#' #====================FGW================
#'
#'
#'
#'alpha <- 0.05
#'
#'p <- c(0.1,0.15,0.1)
#'n <- c(30,15,10)
#'mu <- c(0,0,0)
#'sigma <- c(1,1,1)
#'N <- 1000
#'FGW(n,p,mu,sigma,N)
#'
#'p <- c(0.1,0.15,0.1,0.6)
#'n <- c(30,15,10,50)
#'mu <- c(0,0,0,0)
#'sigma <- c(1,1,1,1)
#'C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#'
#'N <- 1000;
#'FGW(n,p,mu,sigma,N,C2 = C2)
#'
#'
#'
#' #====================FGH===============
#'
#'
#'alpha <- 0.05
#'
#'p <- c(0.1,0.15,0.1)
#'n <- c(30,15,10)
#'mu <- c(0,0,0)
#'sigma <- c(1,1,1)
#'N <- 1000
#'FGH(n,p,mu,sigma,N)

#'p <- c(0.1,0.15,0.1,0.6)
#'n <- c(30,15,10,50)
#'mu <- c(0,0,0,0)
#'sigma <- c(1,1,1,1)
#'C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#'
#'N<-1000;
#'FGH(n,p,mu,sigma,N,C2 = C2)
#'
#'
#'
#' #====================MOVERW=================
#'
#'
#' alpha <- 0.05
#'
#' p <- c(0.1,0.15,0.1)
#' n <- c(30,15,10)
#' mu <- c(0,0,0)
#' sigma <- c(1,1,1)
#' N <- 1000
#'
#' MOVERW(n,p,mu,sigma,N)
#'
#' p <- c(0.1,0.15,0.1,0.6)
#' n <- c(30,15,10,50)
#' mu <- c(0,0,0,0)
#' sigma <- c(1,1,1,1)
#' C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#' N <- 1000
#'
#' MOVERW(n,p,mu,sigma,N,C2 = C2)
#'
#'
#' #====================MOVERH================
#'
#'
#'alpha<-0.05
#'
#'p <- c(0.1,0.15,0.1)
#'n <- c(30,15,10)
#'mu <- c(0,0,0)
#'sigma <- c(1,1,1)
#'N <- 1000
#'MOVERH(n,p,mu,sigma,N)

#'p <- c(0.1,0.15,0.1,0.6)
#'n <- c(30,15,10,50)
#'mu <- c(0,0,0,0)
#'sigma <- c(1,1,1,1)
#'C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#'
#'N <- 1000;
#'MOVERH(n,p,mu,sigma,N,C2 = C2)
#'}
#'
NULL
