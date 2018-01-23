#' The SCIs Based on FGW Methods
#'
#' A method based on generalized pivotal quantity with order statistics(also see \code{\link{FGH}}) to construct the simultaneous confidence intervals for
#' Ratios of Means of Log-normal Populations with excess Zeros.
#'
#' More information about FGW, you can read the paper: Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros.
#'
#' @name FGW
#' @aliases FGW
#' @usage FGW(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)
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
#'FGW(n,p,mu,sigma,N)
#'
#'\dontrun{
#'p <- c(0.1,0.15,0.1,0.6)
#'n <- c(30,15,10,50)
#'mu <- c(0,0,0,0)
#'sigma <- c(1,1,1,1)
#'C2 <- rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
#'
#'N <- 1000;
#'FGW(n,p,mu,sigma,N,C2 = C2)
#'
#'}
#'
#' @export FGW
GPQW0<-function(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05){


  t1<-Sys.time()
  result<-list()

  r = r1 = r2 = matrix(0,dim(C2)[1],N)

  for(i in 1:length(n)){
    assign(paste0('n',i,'0'),rbinom(1,n[i],p[i]))
    assign(paste0('n',i,'1'),n[i]-get(paste0('n',i,'0')))
    assign(paste0('y',i,'bar'),rnorm(1,mu[i],(sigma[i]/sqrt(sqrt(get(paste0('n',i,'1')))))))
    assign(paste0('s',i,'sq'),(sigma[i])^2*rchisq(1,df=(get(paste0('n',i,'1'))-1),ncp=0))

    assign(paste0('Z',i),rnorm(N,0,1))
    assign(paste0('U',i),rchisq(N,df=(get(paste0('n',i,'1'))-1),ncp=0))

    assign(paste0('Z',i,'W'),rnorm(N,0,1))

    #assign(paste0('TPW',i),numeric(N))
    #assign(paste0('TW',i),numeric(N))

  }



  for(k in 1:N){

    tempTW <- list()

    for(i in 1:length(n)){

      assign(paste0('TPW',i,'[',k,']'),
        (get(paste0('n',i,'0')) + 0.5 * (get(paste0('Z',i,'W'))[k])^2)/(n[i] + (get(paste0('Z',i,'W'))[k])^2)-
        (get(paste0('Z',i,'W'))[k] * sqrt(get(paste0('n',i,'0'))*(1 - get(paste0('n',i,'0'))/n[i]) +
          (get(paste0('Z',i,'W'))[k])^2/4))/(n[i] + (get(paste0('Z',i,'W'))[k])^2)

        )
      #TPW1[k]<-(n10+0.5*(Z1W[k])^2)/(n1+(Z1W[k])^2)-(Z1W[k]*sqrt(n10*(1-n10/n1)+(Z1W[k])^2/4))/(n1+(Z1W[k])^2);


      assign(paste0('TW',i,'[',k,']'),log(1-get(paste0('TPW',i,'[',k,']')))+
               (get(paste0('y',i,'bar'))-get(paste0('Z',i))[k]*sqrt(get(paste0('s',1,'sq')))/
               (sqrt(get(paste0('n',i,'1'))*get(paste0('U',i))[k]))+
               get(paste0('s',1,'sq'))/(2*get(paste0('U',i))[k]))
      )

      #TW1[k]<-log(1-TPW1[k])+(y1bar-Z1[k]*sqrt(s1sq)/(sqrt(n11*U1[k]))+s1sq/(2*U1[k]));
      tempTW[[i]] <- get(paste0('TW',i,'[',k,']'))
    }

    r[,k] <- C2 %*% unlist(tempTW)
  }


  for(i in 1: dim(C2)[1]){

    r1[i,] <- sort(r[i,])
    r2[i,] <- rank(r[i,])

  }

  mik1=mak1=numeric(N)

  for(k in 1:N){
    tempr2 <- list()
    for(i in 1: dim(C2)[1]){
      tempr2[[k]] <- r2[i,k]
    }

    mik1[k] <- min(unlist(tempr2))
    mak1[k] <- max(unlist(tempr2))

  }

  mik2<-sort(mik1);
  mak2<-sort(mak1);

  kl<-mik2[N*alpha/2]
  ku<-mak2[N*(1-alpha/2)]

  for( i in 1: dim(C2)[1]){

    assign(paste0('the.',i,'th','lower.','limit'),round(get(paste0('r1'))[i,kl],6))
    assign(paste0('the.',i,'th','upper.','limit'),round(get(paste0('r1'))[i,ku],6))
  }

  result$title1 <- "====================Method: FGW====================="
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
  result$t <- t2-t1;
  result;

}

# result print func
FGW<-function(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)) ,alpha=0.05)
{
  GPQW1<-GPQW0(n,p,mu,sigma,N,C2,alpha);
  for(i in c('title1','title2','interval', 'star',"t"))
  {
    print(GPQW1[[i]])

  }

}











