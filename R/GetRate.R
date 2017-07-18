#' GetOutRates
#'
#'
#' The function calculates two rates - True Positive and False Positive Rates (TPR and FPR) - in order
#' to evaluate the performance
#' of an outlier detection method. Both rates range between zero and one.
#'
#'
#' @param label The numeric (character) vector containing the outlier memberhip, e.g. 1: non-outliers, 0: outliers
#' @param predict The numeric (character) vector consisting of values refering the outlier memberhip predicted by
#' an outlier dection method, e.g. 1: declared as non outliers, 0: identified as outliers
#' @param out a numeric (character) value corresponding to the value indicating outliers, e.g. \code{out=0}
#'
#' @return
#' The function returns list of two values
#' \item{TPR}{The true positive rate  is defined as the ration between the number of
#' correctly idenitified outliers and the number of actual outliers. \code{TPR=1} means that all outliers are correctly identified.}
#' \item{FPR}{The false positive rate is calculated as the number of wrongly identified outliers devided by
#' the total number of actual non-outliers. \code{FPR=0} means that no non-outlier is decleared as outlier.}
#'
#' @examples
#' # generate data
#' d <- SimData(size_grp=c(40,40,40),p_inf=50,
#' p_noise=750,p_out_noise=75)
#' dat <- scale(d$x)
#' lb <- d$lb
#' table(lb) # outliers have zero group memberhsip
#'
#' res <- wrsk(data=dat,k=3,s=6)
#' table(d$lb,res$outclusters)
#'
#' GetOutRates(lb,res$outclusters)
#'
#' @author Sarka Brodinova <sarka.brodinova@tuwien.ac.at>
#'
#' @export
GetOutRates <-function(label,predict,out="out"){
label0 <- rep(FALSE,length(label))
label0[label==out] <- TRUE

predict0 <- rep(FALSE,length(predict))
predict0[predict==out] <- TRUE

tpr <- sum(label0 & predict0)/sum(label0)
fpr <- sum(!label0 & predict0)/sum(!label0)

return(c(TPR=tpr,FPR=fpr))
}


