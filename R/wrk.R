#' wrk
#'
#'
#' This function performs robust (weighted) k-means clustering which is very useful in case of contaminated data. The method
#' aims at detecting both clusters and outliers.
#'
#' @param data A data matrix with n observations and p variables.
#' @param D A distance matrix.
#' @param k The number of clusters.
#' @param max.iter The maximum number of iterations to reach local optimum, the default is 30.
#' @param cutoff A cutoff value to determine outliers. An observation is
#' declared as an outlier if its weight is smaller than or equal to this cutoff, the default is 0.5.
#'
#' @return
#' \item{clusters}{An integer vector with values from 1 to k, indicating a resulting cluster membership.}
#' \item{obsweights}{A numeric vector of observation weights ranging between 0 and 1.}
#' \item{outclusters}{An integer vector with values from 0 to k, containing both cluster membership and identified outliers.
#' 0 corresponds to outlier.}
#' \item{WCSS}{The within-cluster sum of squares of the local optimum.}
#' \item{centers}{The set of cluster centers.}
#'
#' @export
#' @examples
#' # generate data
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#' dat <- scale(d$x)
#'
#' res <- wrk(data=dat,D=dist(dat),k=3)
#' table(d$y,res$clusters)
#' table(d$lb,res$outclusters)
#' plot(res$obsweights,col=d$lb +1)
#'
#' @details
#' \code{wrk} initializes the clustering procedure by ROBIN approach and
#' incorporates a weighting function on each detected clusters during k-means clustering.
#' The weighting function uses LOF in order to assign a weight for each observation. The resulting observation weights reflect
#' a degree of outlyingness and range between 0 and 1. These weights are used to
#' identify outliers as observations with a weight \code{<=cutoff}.
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#'  @references S. Brodinova, P. Filzmoser, T. Ortner, C. Breiteneder, M. Zaharieva. Robust and sparse k-means clustering for
#' high-dimensional data. Submitted for publication, 2017. Available at http://arxiv.org/abs/1709.10012
#'
#' @seealso \code{\link[wrsk]{ROBIN}}
#'
wrk <- function(data,D,k,max.iter=30,cutoff=0.5){
  n <- nrow(data)
  p <- ncol(data)

  # Robin initialization
  init <- ROBIN(D,data,k)
  means <- data[init$centers,] %>% as.matrix()

  go <- TRUE
  iter <- 1
  opt.crit <- +Inf
  oldclass <- rep(1,n)

  while(go){

    # assign points to the closest clusters and calculate their distances
    distToCenter <- matrix(NA,nrow=n,ncol=k)
    for (j in 1:k){
      muj <- means[j,]
      distToCenter[,j] <- apply(scale(data,center=muj,scale=FALSE)^2,1,sum)
    }
    icluster <- apply(distToCenter,1,which.min)
    wcss <- apply(distToCenter,1,min)
    crit <- sum(wcss) # stop criterion

    # calculate observetion weights and cluster centers as weighted means
    A <- GetWeightsAll(n,data,k,icluster,means)
    means <- A$means %>% as.matrix()

    # rearange clustering solution based on the weights
    outcluster <-  icluster
    outcluster[A$obsweights<=cutoff] <- 0


    # local optimal solution
    if( (iter<max.iter) & !identical(oldclass,outcluster) ){
      go <- TRUE
      if(crit<opt.crit){
        #print(table(lb,out.iclass))
        WCSS <- crit
        outCluster <- outcluster
        cluster <-  icluster
        optmeans <- A$means
        optweights <- A$obsweights
      }
      iter <- iter +1
      oldclass <- outcluster
    }else{
      go <- FALSE
      break
    }

  } # end while

  return(list(clusters=cluster,
              obsweights=optweights,
              outclusters=outCluster,
              WCSS=WCSS,
              centers=optmeans))
}




GetWeightsCl <- function(x,lof.k=10){
  # x .... detected cluster
  # lof.k ... the number of nearest neighbors for LOF

  lf <- lof(x,k=lof.k)
  lf.sc <- scale(lf)
  M1 <- median(lf.sc) + 2.5 * mad(lf.sc)
  const1 <- 2
  v1 <- (1 - ((lf.sc - M1)/(const1 - M1))^2)^2
  v1[lf.sc < M1] <- 1
  v1[lf.sc > const1] <- 0

  return(v1)
}



GetWeightsAll <- function(n,data,k,icluster,means){
  v <-numeric(n)
  for(l in 1:k){
    cl <- data[icluster==l,,drop=FALSE] %>% as.matrix()

    if(nrow(cl)>2){
      q <- ifelse(nrow(cl)<=10,nrow(cl)-1,10)
      vo <- GetWeightsCl(x=cl,lof.k=q)
      vo <- as.vector(vo,mode="numeric")
      means[l,] <- apply(vo * cl,2,sum) /sum(vo)
    }else{
      means[l,] <-apply(cl,2,mean)
      vo <- rep(1,nrow(cl))
    }
    v[icluster==l] <- vo
  }
  return(list(obsweights=v,means=means))
}




