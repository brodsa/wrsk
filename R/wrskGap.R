#' wrskGap
#'
#'
#' \code{wrskGap} estimates the optimal sparsity parameter controlling the sparsity in the variable
#' weight vector for \code{wrsk}. The approach uses gap statistic which is additionally adjusted by observation weights.
#' The statistic can be employed to estimate either
#' the sparsity parameter or both the sparsity parameter and the number of clusters in case no information about
#' data is available.
#'
#' @param data A data matrix with n observations and p variables.
#' @param K The number of clusters.
#' @param S A numeric vector containing the values of the sparsity parameter. The values have to be
#' in increasing order, greater than 1 and smaller than \code{sqrt(p)}.
#' @param npermute The number of permutations to calculate gap statistic, the default is 10.
#' @param cores The number of cores for parallel computing, see examples.
#'
#' @return
#' \item{gap}{A numeric vector containing the values of gap statistic for each choice of the sparsity parameter \code{S}.}
#' \item{se}{A numeric vector containing the standard error of log(WCSSS) with respect to the permuted data.}
#' \item{s}{A numeric vector of the investigated values of the sparsity parameter.}
#' \item{resFinal}{A list of objects corresponding to the output of
#' \code{\link[wrsk]{wrsk}}.}
#'
#' @details
#' The estimation of the optimal parameters is based on a data permutation. The gap statistic compares the
#' the adjusted weighted between-cluster sum of squares (WBCSS) obtained by \code{wrsk} with the expected value
#' under a null distribution. The expected value is calculated with respect to WBCSS computed on the clustering solution
#' achieved on the permuted datasets.
#'
#' The optimal sparsity parameter is selected using a 1-standard error rule. If the optimal number of clusters needs to
#' be estimated, the gap statistic is calculated with respect to various choices of \code{S} and {K}.
#' Consequently, for each number of cluster, an optimal sparsity parameter is selected. The optimal number of clusters is chosen based on
#' the values of gap statistic corresponding the optimal sparsity parameter.
#' The largest value of these values indicates the optimal number of clusters.
#'
#'
#' @export
#' @examples
#' \dontrun{
#' # set the number of cores
#' cores <- detectCores() - 2
#'
#' # generate data
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#' dat <- scale(d$x)
#' lb <- d$lb
#'
#' # selection of optimal par. from  K=2,...5 and S=1.1,1.6,..
#' GAP_sk <- rep(-Inf,5)
#' sim_nS <- seq(1.1,sqrt(ncol(dat)),0.5)
#'
#' # --- K=2 ------
#' k2 <- wrskGap(data=dat,K=2,S=sim_nS,npermut=10,cores)
#' id.max.k2 <- which.max(k2$gap)
#' id.opt.k2 <- which(k2$gap[1:id.max.k2]> (k2$gap[id.max.k2] - k2$se[id.max.k2]))[1]
#' GAP_sk[2] <- k2$gap[id.opt.k2]
#'
#' # ---- K=3 ------------
#' k3 <- wrskGap(data=dat,K=3,S=sim_nS,npermut=10,cores)
#' id.max.k3 <- which.max(k3$gap)
#' id.opt.k3 <- which(k3$gap[1:id.max.k3]> (k3$gap[id.max.k3] - k3$se[id.max.k3]))[1]
#' GAP_sk[3] <- k3$gap[id.opt.k3]
#'
#' # ---- K=4 ------------
#' k4 <- wrskGap(data=dat,K=4,S=sim_nS,npermut=10,cores)
#' id.max <- which.max(k4$gap)
#' id.opt.k4 <- which(k4$gap[1:id.max]> (k4$gap[id.max] - k4$se[id.max]))[1]
#' GAP_sk[4] <- k4$gap[id.opt.k4]
#'
#' # --- K=5 -----------
#' k5 <- wrskGap(data=dat,K=5,S=sim_nS,npermut=10,cores)
#' id.max <- which.max(k5$gap)
#' id.opt.k5 <- which(k5$gap[1:id.max]> (k5$gap[id.max] - k5$se[id.max]))[1]
#' GAP_sk[5] <- k5$gap[id.opt.k5]
#'
#'
#' # THE OPTIMAL CLUSTERING SOLUTION
#' id.opt.k <- which.max(GAP_sk)
#' opt.k.solution <- list(0,k2,k3,k4,k5,k6)[[id.opt.k]]
#' id.opt.s <- c(0,id.opt.k2,id.opt.k3,id.opt.k4,id.opt.k5,id.opt.k6)[id.opt.k]
#' opt.res <- opt.k.solution$resFinal[[id.opt.s]]
#'
#' obsweights <- opt.res$obsweights
#' varweights <- opt.res$varweights
#' clusters <- opt.res$clusters
#' outclusters <- opt.res$outclusters
#'
#' table(d$lb,outclusters)
#' plot(varweights)
#' plot(obsweights,col=d$lb+1)
#' }
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#' @references S. Brodinova, P. Filzmoser, T. Ortner, C. Breiteneder, M. Zaharieva. Robust and sparse k-means clustering for
#' high-dimensional data, preparing for submission, 2017.
#' 
#' @references D. M. Witten and R. Tibshirani. A framework for feature selection in clustering.
#' Journal of the American Statistical Association, 105(490) 713-726, 2010.
#'
#' @seealso \code{\link[wrsk]{wrsk}}
#'
#' @importFrom parallel detectCores
#'
wrskGap <- function(data,K,S,npermute=10,cores){
  rownames(data) <- 1:nrow(data)
  colnames(data) <- 1:ncol(data)

  gap_s <- se_s <- logBsa <- logBs <- c()
  B_s <- nonZerW <-  c()
  resFinal <- list()
  for(j in 1:length(S)){
    cat(S[j])
    #ress <- RobWeightSparseKmeans(dat=data,k=K,s=S[j],weights_par=weights_par_our)
    ress <- tryCatch(wrsk(dat=data,k=K,s=S[j]),
             warning=function(w){print(paste("warnings"))},
             error=function(e){print(paste("one cluster disappeared - the procedure restared at s = ",S[j]))})

    if(ress[1] == paste("one cluster disappeared - the procedure restared at s = ",S[j])){
      ress <- tryCatch(wrsk(dat=data,k=K,s=S[j]),
                       warning=function(w){print(paste("warnings"))},
                       error=function(e){print(paste("one cluster disappeared - the procedure stopped at s = ",S[j]))})

      if(ress[1] == paste("one cluster disappeared - the procedure stopped at s = ",S[j])){
        return(list(gap=gap_s,se=se_s,
                    nonZerW=nonZerW,s=S[1:j],resFinal=resFinal))
      }
    }
    resFinal[[j]] <- ress
    B_s <- ress$WBCSS
    nonZerW <- c(nonZerW,length(which(ress$varweights!=0)))

    B_sa <- numeric(npermute)
    res.per <-mclapply(1:npermute,
                       function(x,dat,k,s,weights_par){
                         set.seed(x)
                         dat.per <- apply(dat,2,function(z){sample(z,length(z))})
                         r <- wrsk(dat.per,k,s)
                       },dat=data,k=K,s=S[j],mc.cores=cores)

    B_sa <- numeric(npermute)
    for(a in 1:npermute){
      B_sa[a] <- res.per[[a]]$WBCSS
    }

    gap_s <- c(gap_s,log(B_s)- mean(log(B_sa)))
    se_s <- c(se_s,1*(  sqrt(1+1/length(B_sa))*sqrt(var(log(B_sa)))*sqrt((length(B_sa)-1)/length(B_sa)))  )

    if(nonZerW[j]==ncol(data)| S[j]==tail(S,1)){
      break
      return(list(gap=gap_s,se=se_s,
                  nonZerW=nonZerW,s=S[1:j],resFinal=resFinal))
    }
  }
  return(list(gap=gap_s,se=se_s,
              nonZerW=nonZerW,s=S[1:j],resFinal=resFinal))

}

