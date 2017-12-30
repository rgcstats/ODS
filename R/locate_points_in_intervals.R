#' Identifies which intervals some values lie within.
#' @param y numeric vector
#' @param cutoffs A set of finite cutoffs defining the intervals.
#'  If there are H cutoffs there are H+1 intervals
#' @return a list with elements:
#' $interval.number: a vector of n integers (where n=length(y)) defining which
#' interval contains #' each element of y. A value of 1 corresponds to the
#' interval [-Inf, cutoff[1]), 2 corresponds to the interval
#'  [cutoff[1],cutoff[2]), ..., H+1 corresponds to [cutoff[H],Inf).
#' $lower.cutoffs: the lower bounds of the H+1 intervals
#' $upper.cutoffs: the upper bounds of the H+1 intervals
#' @examples
#' test <- locate.points.in.intervals(y=c(-100,-5:5),cutoffs=c(-10,-2,2,4,10))
#' cbind(y=c(-100,-5:5),interval=test$interval.number,
#'       lower=test$lower.cutoffs,
#'       upper=test$upper.cutoffs)
#' @export
locate.points.in.intervals <- function(y,cutoffs){
  # cutoffs must be a vector of finite interval cutoffs
  #  if there are H intervals there should be (H-1) cutoffs
  cutoffs <- c(-Inf,sort(cutoffs),Inf)
  cutoff.le <- outer(cutoffs[-length(cutoffs)],y,"<=")
  cutoff.gt <- outer(cutoffs[-1],y,">")
  lower.cutoff.indexes <- apply(cutoff.le*cutoff.gt,2,which.max)
  return(list(interval.number=lower.cutoff.indexes,
              lower.cutoffs=cutoffs[lower.cutoff.indexes],
              upper.cutoffs=cutoffs[lower.cutoff.indexes+1]))
}
