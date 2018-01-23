#' The SCIs Based on FGH Methods
#'
#' A method based on generalized pivotal quantity with order statistics(also see \code{\link{FGW}}) to construct the simultaneous confidence intervals for
#' Ratios of Means of Log-normal Populations with excess Zeros.
#'
#' More information about FGH, you can read the paper: Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros.
#'
#' @name FGH
#' @aliases FGH
#'
#' @usage FGH(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)
#'
#' @param n The sample size of the mixture distributions,must be an integer vector.
#' @param p The zero probability of the mixture distribution,it has the same length to the \strong{n} params.
#' @param mu The mean of the non-zero samples,which after log-transformation.
#' @param sigma The variance of the non-zero samples,which after log-transformation.
#' @param N The number of independent generated data sets.
#' @param C2 Matrix C,You can refer to the paper of Xu et al. for specific forms.
#' @param alpha The confidence level,it always set \emph{alpha=0.5}
#'
#' @return The method will return the Simultaneous Confidence Intervals(SCIs) and the time consuming
#'
#' @author Jing Xu, Xinmin Li, Hua Liang
#'
#' @examples
#'
#'
#'alpha <- 0.05
#'
#'p <- c(0.1,0.15,0.1)
#'n <- c(30,15,50)
#'mu <- c(0,0,0)
#'sigma <- c(1,1,1)
#'N <- 500
#'FGH(n,p,mu,sigma,N)
#'
#'\dontrun{
#'p <- c(0.1,0.15,0.1,0.6)
#'n <- c(30,15,10,50)
#'mu <- c(0,0,0,0)
#'sigma <- c(1,1,1,1)
#'C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#'
#'N <- 1000;
#'FGH(n,p,mu,sigma,N,C2 = C2)
#'}
#'
#'
#' @export FGH

rndMixture<-function(N0,N1,p,N){

  u<-numeric(N);
  u<-runif(N,0,1);
  r<-rep(0,N);
  for(i in 1:N)
  {
    if(u[i]<p)
    {r[i]<-rbeta(1,N0+1,N1);}
    else{r[i]<-rbeta(1,N0,N1+1);}
  }
  r
}

GPQH0<-function(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05){

  t1<-Sys.time()
  result<-list()

  r = r1 = r2 = matrix(0,dim(C2)[1],N)

  for(i in seq_along(n)){

    assign(paste0('n',i,'0'),rbinom(1,n[i],p[i]))
    assign(paste0('n',i,'1'),n[i]-get(paste0('n',i,'0')))
    assign(paste0('y',i,'bar'),rnorm(1,mu[i],sigma[i]/sqrt(get(paste0('n',i,'1')))))
    assign(paste0('s',i,'sq'),sigma[i]^2*rchisq(1,df=(get(paste0('n',i,1))-1),ncp=0))

    assign(paste0('Z',i),rnorm(N,0,1))
    assign(paste0('U',i),rchisq(N,df=(get(paste0('n',i,'1'))-1),ncp=0))

    assign(paste0('TPH',i),rndMixture(get(paste0('n',i,'0')),get(paste0('n',i,'1')),0.5,N))

    #assign(paste0('TPW',i),numeric(N))
    #assign(paste0('TW',i),numeric(N))
  }


  for(k in 1:N){

    tempTW <- list()

    for(i in seq_along(n)){
      assign(paste0('TW',i,'[',k,']'),log(1-get(paste0('TPH',i))[k])+
        (get(paste0('y',i,'bar')) - get(paste0('Z',i))[k]*sqrt(get(paste0('s',i,'sq')))/
          (sqrt(get(paste0('n',i,'1'))*get(paste0('U',i))[k]))+
          get(paste0('s',i,'sq'))/(2*get(paste0('U',i))[k]))
        )
      #TW1[k]<-log(1-TPH1[k])+(y1bar-Z1[k]*sqrt(s1sq)/(sqrt(n11*U1[k]))+s1sq/(2*U1[k]));

      tempTW[[i]] <- get(paste0('TW',i,'[',k,']'))

      }

    r[,k]<-C2%*%unlist(tempTW)

  }

  for(i in 1:dim(C2)[1]){
    r1[i,] <- sort(r[i,])
    r2[i,] <- rank(r[i,])
  }


  mik1=mak1=numeric(N)

  for(k in 1:N){
    tempr2 <- list()
    for(i in 1:dim(C2)[1]){

      tempr2[[i]] <- r2[i,k]
    }
    mik1[k] <- min(unlist(tempr2))
    mak1[k] <- max(unlist(tempr2))
  }

  mik2<-sort(mik1);
  mak2<-sort(mak1);

  kl<-mik2[N*alpha/2]
  ku<-mak2[N*(1-alpha/2)]

  for( i in 1: dim(C2)[1]){

    assign(paste0('the.',i,'th','lower.','limit'),round(r1[i,kl],6))
    assign(paste0('the.',i,'th','upper.','limit'),round(r1[i,ku],6))
  }


  # result$the.first.lower.limit<-r1[1,kl];
  # result$the.first.upper.limit<-r1[1,ku];
  #
  # result$the.second.lower.limit<-r1[2,kl];
  # result$the.second.upper.limit<-r1[2,ku];
  #
  # result$the.third.lower.limit<-r1[3,kl];
  # result$the.third.upper.limit<-r1[3,ku];

  result$title1 <- "====================Method: FGH====================="
  result$title2 <- "The Simultaneous Confidence Intervals are:          "

  tempinterval <- list()
  for(i in 1: dim(C2)[1]){
    tempinterval[[i]] <- paste0('[',get(paste0('the.',i,'th','lower.','limit')),',',get(paste0('the.',i,'th','upper.','limit')),']')

  }

  interval0 <- do.call(rbind,tempinterval)
  interval1 <- data.frame(interval0)
  names(interval1) <- c('[LCL,UCL]')

  result$interval <- interval1

  t2 <- Sys.time()
  result$star <- '**********************Time**************************'
  result$t <- t2-t1
  result

}


FGH<-function(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)
{
  GPQH1<-GPQH0(n,p,mu,sigma,N,C2,alpha);
  for(i in c('title1','title2','interval', 'star',"t"))
  {
    print(GPQH1[[i]])

  }

}




