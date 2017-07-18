#' SimData
#'
#' @description
#' This function generates a synthetic dataset for the purpose of data clustering
#' in case of both outliers and noise variables are present.
#'
#' @param size_grp A numeric vector containing the group sizes to be generated. The length of the vector
#' corresponds to the number of groups to be generated.
#' @param p_inf The number of informative variables describing the generated group structure.
#' @param p_noise The number of noise variables which do not contribute to the group separation. If
#' not specified, no noise variables are generated.
#' @param p_out_inf The number of informative variables in which the observations are contaminated,
#' i.e. replaced by outliers either scatter outliers or uniformly distributed outliers, see \code{scatter_out}.
#' @param pct_out The proportion of observations to be contaminated in the informative variables, default is 0.10.
#' @param scatter_out If \code{TRUE}, scattered outliers are generated with the characteristics
#' specified in \code{s_out_range}, otherwise uniformly distributed outliers are produced with the
#' specification defined in \code{unif_out_range}.
#' @param p_out_noise The number of noise variables in which the contamination is conducted,
#' see \code{unif_out_range}.
#' @param noise_pct_out The proportion of observations to be contaminated in noise variables by replacing
#' them with uniformly distributed outliers. The contaminated observations differ
#' from those contaminated in informative variables, default is 0.10.
#' @param unif_out_range Optional argument. You can change the interval of an uniform distribution
#' to randomly generated outliers in \code{[min1,max1] or [min2,max2]}. The specification hase to be in the list,
#' default is \code{list(min1=-12,max1=-6,min2=6,max2=12)}.
#' @param mu_grp_range Optional argument, see references. Default is \code{list(min1=-6,max1=-3,min2=3, max2=6)}.
#' @param s_out_range Optional argument, see references. Default is \code{list(min=3,max=10)}.
#' @param rho_grp_range Optional argument, see references. Default is \code{list(min=0.1,max=0.9)}.
#'
#' @return
#' \item{x}{A data matrix of a synthetic dataset.}
#' \item{y}{An integer vector corresponding to a  group membership before contamination.}
#' \item{lb}{An integer vector with group labels and outlier labels denoted by 0.}
#' \item{lbout}{An integer vector with group lables and labels for outlier in informative variables (0)
#' and noise variables (the number of groups+1).}
#'
#' @export
#'
#' @details
#' Groups are generated in the first \code{p_inf} informative variables with various characteristics
#' following Gaussian models. The groups have different mean vectors and covariance matrix which is additionally
#' randomly rotated. If uninformative variables are required, \code{p_noise} noise variables
#' are generated following uniform distributions and added to an informative part.
#' Two types of outliers - scattered and uniformly distributed- are considered to contaminated data.
#' The outiers can be placed either in the informative or uninformative part.
#'
#' @examples
#' # Generate 3 groups of equal sizes in the first 50 variables with 10% of
#' # scatter outliers in all 50 informative variables, and
#' # 10% of uniformly distributed outliers in 75 noise variables.
#'
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#'
#' # group membership with outliers in 0 group
#' table(d$lb)
#'
#' # scatter outliers in 0 group and uniformly distributed outliers in 4 group
#' table(d$lbout)
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#' @references S. Brodinova, P. Filzmoser, T. Ortner, C. Breiteneder, M. Zaharieva. Robust and sparse k-means clustering for
#' high-dimensional data, preparing for submission, 2017.
#'
#' @importFrom mixAK rRotationMatrix
#' @importFrom MASS mvrnorm



SimData <- function(size_grp,p_inf,p_noise=NULL,p_out_inf=NULL,
                    pct_out=0.10,scatter_out=TRUE,p_out_noise=NULL,noise_pct_out=0.10,
                    unif_out_range=NULL,mu_grp_range=NULL,s_out_range=NULL,
                    rho_grp_range=NULL){
  if(is.null(unif_out_range)){
    unif_out_range <- list(min1=-12,max1=-6,
                           min2=6,max2=12)
  }
  if(is.null(mu_grp_range)){
    mu_grp_range <- list(min1=-6,max1=-3,
                         min2=3, max2=6)
  }
  if(is.null(s_out_range)){
    s_out_range <- list(min=3,max=10)
  }
  if(is.null(rho_grp_range)){
    rho_grp_range <-list(min=0.1,max=0.9)
  }

  GenSynData<- function(size_grp,p_inf,p_noise,p_out_inf,
                           n_grp,s_grp,mu_grp,
                           n_out,s_out,scatter_out,
                        unif_out_range,
                        mu_grp_range,
                        s_out_range,rho_grp_range){



    # caluculate center of groups
    if(p_inf==2){
      MU_grp <- matrix(0,nrow=p_inf+1, ncol=n_grp)
      diag(MU_grp) <- mu_grp
      MU_grp <- MU_grp[1:2,]
    }else{
      n_p <- floor(p_inf/n_grp)
      MU_grp <-  matrix(nrow=0, ncol=n_grp)
      for(i in 1:n_p){
        MU_grp <- rbind(matrix(0,nrow=n_grp, ncol=n_grp),MU_grp)
        diag(MU_grp) <- mu_grp
      }
      if(nrow(MU_grp)!=p_inf){
        for(i in 1:(p_inf-nrow(MU_grp))){
          mi <- rep(0,n_grp)
          mi[i] <- mu_grp[i]
          MU_grp <- rbind(MU_grp,mi)
        }
      }
    }

    # the outliers has the same location as groups (only in case of scatter outliers)
    mu_out <- MU_grp
    n <- sum(size_grp)
    p <- p_inf + p_noise
    lb1 <- vector(length=0)
    y <- vector(length=0)
    x <- matrix(nrow=0, ncol=p)


    for(i in 1:n_grp){
      lb <- c(rep(i,size_grp[i]))
      y <- c(y,rep(i,size_grp[i]))
      S_grp <- s_grp[[i]]

      X <- mvrnorm(size_grp[i],mu=MU_grp[,i], Sigma=S_grp)
      # noise variables
      if(!is.null(p_noise)){
        X <- cbind(X, matrix(rnorm(p_noise*size_grp[i], mean=0), ncol=p_noise))
      }
      # generate outliers in informative variables
      if(n_out[i]!=0){
        if(is.null(p_out_inf)){
          p_out <- 1:p_inf
        }else{
          p_out <- 1:p_out_inf
        }
        S_out <- matrix(0, ncol=length(p_out), nrow=length(p_out))
        diag(S_out) <- s_out[i]

        # scatter outliers
        if(scatter_out){
          X[1:n_out[i],p_out] <-mvrnorm(n_out[i],mu=mu_out[p_out,i], Sigma=S_out)
          lb[1:n_out[i]] <-0

        # non-scatter outliers (outliers with differnt location than groups)
        }else{
          X[1:n_out[i],p_out] <-matrix(sample(c(runif(n_out[i]*length(p_out),
                                                      min=unif_out_range$min1,max=unif_out_range$max1),
                                                runif(n_out[i]*length(p_out),
                                                      min=unif_out_range$min2,max=unif_out_range$max2)),
                                              n_out[i]*length(p_out)),
                                         nrow=n_out[i],ncol=length(p_out))
          lb[1:n_out[i]] <-0
        }

      }
      x <- rbind(x,X)
      lb1 <- c(lb1,lb)
    }
    colnames(x) <- paste0("x", rep(1:p))
    return(list(x=x,y=y,lb=lb1))
    # x ... a data set
    # y ... group membership before contamination
    # lb ... group and outlier membership
  }

  # groups
  n_grp <-  length(size_grp)
  mu_grp <- rep(sample(c(runif(n_grp,min=mu_grp_range$min1,max=mu_grp_range$max1),
                         runif(n_grp,min=mu_grp_range$min2,max=mu_grp_range$max2)),1),n_grp)
  ss_grp <- rep(1,n_grp)
  rrho_grp <- runif(n_grp,min=rho_grp_range$min,max=rho_grp_range$max)

  # group covariance structure
  s_grp <- list()
  for( j in 1: n_grp){
    S_g <- matrix(rrho_grp[j], ncol=p_inf, nrow=p_inf)
    diag(S_g) <- ss_grp[j]
    R <- rRotationMatrix(n=1, dim=p_inf)

    s_grp[[j]] <- R %*% S_g %*% t(R)
  }

  # scatter outliers
  n_out <- c(round(size_grp*pct_out))
  rho_out <- 0
  s_out <-runif(n_grp,min=s_out_range$min,max=s_out_range$max)


  ds <- GenSynData(size_grp,p_inf,p_noise,p_out_inf,
                    n_grp,s_grp,mu_grp=mu_grp,
                    n_out,s_out,scatter_out,unif_out_range=unif_out_range,
                   mu_grp_range=mu_grp_range,
                   s_out_range=s_out_range,rho_grp_range=rho_grp_range)

  # outliers in non informative variables
  if(!is.null(p_out_noise)){
    out_p <- sample((p_inf+1):(p_noise+p_inf),p_out_noise)

    a_n_out <- c(round(size_grp*noise_pct_out))
    oo <- c(c(a_n_out,0)+c(0,cumsum(size_grp))+1)[1:n_grp]

    out_o <- apply(cbind(oo,cumsum(size_grp),a_n_out),1,
                   function(x){
                     sample(c(x[1]:x[2]),x[3])
                   })
    out_o <- unlist(c(out_o))
    ds$x[out_o,out_p] <- matrix(sample(c(runif(length(out_o)*length(out_p),
                                               min=unif_out_range$min1,max=unif_out_range$max1),
                                         runif(length(out_p)*length(out_o),
                                               min=unif_out_range$min2,max=unif_out_range$max2)),
                                        length(out_p)*length(out_o)),
                                 nrow=length(out_o),ncol=length(out_p))

    ds$lb[out_o] <- 0

    lb_out <- ds$lb
    lb_out[out_o] <- n_grp+1
  }
  else{
    lb_out <- ds$lb
  }
  return(list(x=ds$x,y=ds$y,lb=ds$lb,lbout=lb_out))

}

