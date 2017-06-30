#' wrsk
#'
#'
#' This function performs robust (weighted) and sparse k-means clustering for high-dimensional data (Brodinova et al (2017)).
#' For the given number of clusters \code{k} and the sparsity parameter \code{s}, the algorithm
#' detects clusters, outliers, and informative variables simultaneously.
#'
#'
#' @param data A data matrix with n observations and p variables.
#' @param k The number of clusters.
#' @param s The sparsity parameter which penalizes the L1 norm of variable weights, i.e. lasso type penalty.
#' The value should be larger than 1 and smaller than \code{sqrt(p)}.
#' @param iteration The maximum number of iterations allowed.
#' @param cutoff A cutoff value to determine outliers. An observation is
#' declared as an outlier if its weight is smaller than or equal to this cutoff, the default is 0.5.
#'
#' @return
#' \item{clusters}{An integer vector with values from 1 to k, indicating a resulting cluster membership.}
#' \item{obsweights}{A numeric vector of observation weights ranging between 0 and 1.}
#' \item{outclusters}{An integer vector with values from 0 to k, containing both cluster membership and identified outliers.
#' 0 corresponds to outlier.}
#' \item{varweights}{A numeric vector of variable weights reflecting the contribution of variables to a cluster separation.
#' A high weight suggests that a variable is informative.}
#' \item{WBCSS}{The weighted-between cluster sum of squares for the local optimum. The value
#' is calculated with respect to the final variable weights and adjusted by the final observation weights.}
#' \item{centers}{The set of final cluster centers.}
#'
#' @export
#'
#' @examples
#' # generate data
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#' dat <- scale(d$x)
#'
#' res <- wrsk(data=dat,k=3,s=6)
#' table(d$lb,res$outclusters)
#' plot(res$varweights)
#'
#' @details
#' The method is a three-step iterative procedure. First, a weighting function is employed
#' during sparse k-means clustering with ROBIN initialization. Then, the variable weights
#' from sparse k-means are updated for the given sparsity parameter. These two steps are
#' repeated until the variable weights stabilize. Finally, both clusters and outliers are detected.
#' The approach is a robust version of sparse k-means (Witten and Tibshirani, 2010) and an alternative of
#' robust (trimmed) and sparse k-means (Kondo et al, 2016).
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#' @references S. Brodinova, P. Filzmoser, T. Ortner, C. Breiteneder, M. Zaharieva. Robust and sparse k-means clustering for
#' high-dimensional data, 2017.
#' @references D. M. Witten and R. Tibshirani. A framework for feature selection in clustering.
#' Journal of the American Statistical Association, 105(490), 713-726, 2010.
#' @references Y. Kondo, M. Salibian-Barrera, R.H. Zamar. RSKC: An R Package for a
#' Robust and Sparse K-Means Clustering Algorithm., Journal of Statistical Software, 72(5), 1-26, 2016.
#'
#' @seealso \code{\link[wrsk]{wrk}}, \code{\link[sparcl]{KmeansSparseCluster}}, \code{\link[RSKC]{RSKC}}
#'
wrsk <- function(data,k,s,iteration=15,cutoff=0.5){
  p <- ncol(data)
  n <- nrow(data)

  oldss<-(-Inf)
  Wp<-rep(1/sqrt(p),p) # initial variable weights

  for(r in 1:iteration){
    ### Step A: Calculation of observation weights for given weights
    ## iteration > 1
    w.data<-t(t(data[,Wp!=0,drop=FALSE])*sqrt(Wp[Wp!=0]))
    if(r>1){
      means <- matrix(NA,nrow=k,ncol=ncol(w.data))
      for (l in 1:k){
        cl <- w.data[icluster==l,,drop=FALSE] %>% as.matrix()
        means[l,] <- apply(Wo[icluster==l] * cl,2,sum) /sum(Wo[icluster==l])
      }
    }

    wD <- dist(w.data)
    # A1 obs weights on weighted data with Wj (variables)
    A1<- wrk(data=w.data,D=wD,k=k,cutoff=cutoff)
    icluster <- A1$clusters
    crit.A1 <- A1$WCSS

    ## A2 obs weights on unweighted data
    A2 <- GetWeightsAll(n=n,data=data,k=k,icluster=icluster,
                        means=matrix(NA,ncol=ncol(data),nrow=k))

    # combination of weights for observation from A1 and A2
    Wo <-apply(cbind(A2$obsweights,A1$obsweights),1,min)

    ### STEP B: update variable weights
    B <- step.b(s=s,data=data,k=k,icluster=icluster,Wp=Wp,Wo=Wo)
    Wp <- B$Wp
    newss <- B$WBSS

    # when local optimal founded
    if(((newss - oldss)/newss) < 1e-8 | r==iteration){
      # last iteration
      w.data<-t(t(data[,Wp!=0,drop=FALSE])*sqrt(Wp[Wp!=0])) # reduce the dimention so that the alg is more efficient
      means <- matrix(NA,nrow=k,ncol=ncol(w.data))
      for (l in 1:k){
        cl <- w.data[icluster==l,,drop=FALSE] %>% as.matrix()
        means[l,] <- apply(Wo[icluster==l] * cl,2,sum) /sum(Wo[icluster==l])
      }

      wD <- dist(w.data)
      # A1 weighted data with Wj (variables)
      A1<- wrk(data=w.data,D=wD,k=k,cutoff=cutoff)
      icluster <- A1$clusters

      # A2 unweighted data
      A2 <- GetWeightsAll(n=n,data=data,k=k,icluster=icluster,
                          means=matrix(NA,ncol=ncol(data),nrow=k))

      # combination of weights for observation from A1 and A2
      Wo <-apply(cbind(A2$obsweights,A1$obsweights),1,min)

      # identifying outliers based final weights for observation and cut of value equal to 0.5
      foutclusters <- icluster      # cluster (>0) and outlier (0) membership
      foutclusters[Wo<=cutoff] <- 0
      fobsweights <- Wo     # weights for observations
      fclusters <- icluster   # cluster membership only
      fvarweights <- Wp      # weights for variables
      fbcss <- B$WBSS
      #print(table(lb,fmemberhip))
      break
    }else{
      oldss <- newss
      r <- r +1
    }
  } # end for loop with repsect to r-iteration

  return(list(clusters=fclusters,
              outclusters=foutclusters,
              obsweights=fobsweights,
              varweights=fvarweights,
              WBCSS=fbcss))
}


step.b<-
  function(s,data,k,icluster,Wp,Wo){ # d is indexed data
    a<-step.b.BSS_js(data=data,icluster=icluster,k=k,Wp=Wp,Wo=Wo)
    a<-pmax(a,0)
    de<-BinarySearch.delta(tss.wss=a,s=s)
    SS<-pmax(a-de,0) # p by 1
    Sn2<-norm2(SS)
    Wp<-SS/Sn2
    WBSS<-sum(Wp*a)
    return(list(Wp=Wp,WBSS=WBSS))
  }


step.b.BSS_js<-function(data,icluster,k,Wp,Wo){
  a<-matrix(NA,nrow=k,ncol=ncol(data))
  sumWp<-sum(Wp)
  Tmu<-apply(Wo * data,2,sum) /sum(Wo)
  Tterm<- Wo * scale(data,center=Tmu,scale=FALSE)^2
  for (i in 1:k){
    if (sum(icluster==i)==0){
      next;
    }else{
      datak<-data[icluster==i,,drop=FALSE] #nk by p
      Wof <- Wo[icluster==i]
      Wmu<-apply(Wof * datak,2,sum) /sum(Wof) # p by 1
      Wterm<-Wof * scale(datak ,center=Wmu,scale=FALSE)^2 #nk by p

      Ttermk<-Tterm[icluster==i,,drop=FALSE] # nk by p
      T_W<-Ttermk-Wterm # nk by p

      adjust <- as.vector(sumWp / ( (!is.na(Wterm)) %*% Wp ))
      adjustTerm<-T_W*adjust #nk by p
      a[i,]<-colSums(adjustTerm,na.rm=TRUE) #p by 1
    }
  }
  a<-colSums(a,na.rm=TRUE)
  return(a)
}



BinarySearch.delta <- function(tss.wss,s){
  # L1=l1bound
  # tss.wss= tssperfeatre - wssperfeature
  if(norm2(tss.wss)==0 || sum(abs(tss.wss/norm2(tss.wss)))<=s) return(0)
  lam1 <- 0
  lam2 <- max(abs(tss.wss))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(tss.wss,(lam1+lam2)/2)
    if(sum(abs(su/norm2(su)))<s){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

norm1 <-function(y){sum(abs(y))}

norm2 <-function(y){sqrt(sum(y^2))}
