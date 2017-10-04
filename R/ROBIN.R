#' ROBIN (ROBust INitialization)
#'
#' @description
#' \code{ROBIN} searches for k initial cluster seeds for k-means-based clustering methods.
#'
#' @param D A distance matrix calculated on \code{data}.
#' @param data A data matrix with n observations and p variables.
#' @param k The number of cluster centers to find.
#' @param mp The number of the nearest neighbors to find dense regions by LOF, the default is 10.
#' @param critRobin The cutoff value for LOF to determine the observations in the dense regions,
#' the default is 1.05.
#'
#'
#' @return
#' \item{centers}{A numeric vector of \code{k} initial cluster centers corresponding to the k indices of observations.}
#' \item{lof}{A real vector of local outlier factor values.}
#'
#' @export
#'
#' @details
#' The centers are the observations located in the most dense region and far away from each other at the same time.
#' In order to find the observations in the highly dense region, ROBIN uses LOF
#' (Local Outlier Factor, Breunig et al (2000)), see more details \code{\link[dbscan]{lof}}.
#'
#' @examples
#' # generate data
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#' dat <- scale(d$x)
#'
#' res <- ROBIN(D=dist(dat),data=dat,k=3)
#' km <- kmeans(x=dat,centers=dat[res$centers,])
#' table(d$y,km$cluster)
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#' @references Hasan AM, et al. Robust partitional clustering by
#' outlier and density insensitive seeding. Pattern Recognition Letters, 30(11), 994-1002, 2009.
#'
#' @seealso \code{\link[dbscan]{lof}}, \code{\link[wrsk]{wrk}}
#'
#' @importFrom dbscan lof
#' @importFrom dplyr %>%
#'
#'
ROBIN <- function(D,data,k,mp=10,critRobin=1.05){
  lf <-lof(D,k=mp)
  n <- nrow(data)
  r <-sample(n,1)
  id.means <- numeric(k)
  m <-1
  D <- as.matrix(D)
  while(m<=k){
    # find the observations that are far away from each other - the potential cluster seeds
    if(m<=2){
      sort.points <- sort.int(D[r,],decreasing = TRUE) %>%
        names() %>% as.numeric()
    }else{
      sort.points <- sort.int(apply(D[id.means[1:m],],2,min),decreasing=TRUE) %>%
        names() %>% as.numeric()
    }
    lf.sort.points <- lf[sort.points]
    id <- which(c(lf.sort.points<critRobin))[1]
    r <- sort.points[id]
    id.means[m] <- r

    m <- m+1
  }
  return(list(centers=id.means,lof=lf))
}


