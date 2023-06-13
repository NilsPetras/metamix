#' t Fit
#'
#' Test the fit of the model against the observed distribution of t-values.
#'
#' @param t observed t-values
#' @param n1 observed sample size of group 1
#' @param n2 observed sample size of group 2
#' @param alpha assumed alpha of observed t-tests
#' @param p ml estimate of the probability of selective publishing
#' @param d ml estimate of the true effect size
#' @param nrep number of bootstrap repetitions (total number will be nrep *
#'   length(n1))
#' @param TwoSided assumption: were the observed t-tests two-sided?
#' @param SigSuppress assumption: were significant results suppressed?
#'
#' @return list including the ks-test output, a matrix of the bootstrapped
#'   t-values, and a matrix of the bootstrapped decision to publish / not
#'   publish
tFit <- function(
    t, n1, n2, alpha, p, d, nrep, TwoSided, SigSuppress
) {
  # if n2 = NULL one sample t test
  
  # the non-mixture model is essentially a mixture model with p_sp = 1
  # this is already passed on to this function by the metamix function
  
  tdist <- matrix(nrow = nrep, ncol = length(n1))
  
  df <- n1 + n2 - 2
  
    for (j in 1:length(n1)) {
      tdist[ ,j] <- stats::rt(
        nrep,
        df[j],
        d * sqrt(n1[j] * n2[j] / (n1[j] + n2[j]))) # see Faul et al. (2007), G-Power paper, Table 3
    }
  
  published <- matrix(
    stats::runif(nrep * length(n1)),
    nrow = nrep,
    ncol = length(n1))
  
  # decision tree with 4 nodes based on:
  # 1) are significant results suppressed? (no/yes)
  # 2) are tests conducted two-sided? (no/yes)
  if (!SigSuppress) {
    if (!TwoSided) {
      tcrit <- stats::qt(1 - alpha, df)
        for (j in 1:length(n1)) {
          published[tdist[,j] > tcrit[j],j] <- TRUE
          published[tdist[,j] < tcrit[j],j] <- published[tdist[,j] < tcrit[j],j] > p
        }
    } else {
      tcrit <- stats::qt(1 - alpha / 2, df)
        for (j in 1:length(n1)) {
          published[abs(tdist[,j]) > tcrit[j],j] <- TRUE
          published[abs(tdist[,j]) < tcrit[j],j] <- published[abs(tdist[,j]) < tcrit[j],j] > p
        } 
    }
  } else {
    if (!TwoSided) {
      tcrit <- stats::qt(1 - alpha, df)
        for (j in 1:length(n1)) {
          published[tdist[,j] < tcrit[j],j] <- TRUE
          published[tdist[,j] > tcrit[j],j] <- published[tdist[,j] > tcrit[j],j] > p
        }
    } else {
      tcrit <- stats::qt(1 - alpha / 2, df)
        for (j in 1:length(n1)) {
          published[abs(tdist[,j]) < tcrit[j],j] <- TRUE
          published[abs(tdist[,j]) > tcrit[j],j] <- published[abs(tdist[,j]) > tcrit[j],j] > p
        }
    }
  }
  
  
  ks <- stats::ks.test(t, c(tdist[as.logical(published)]))
  ks$data.name <- "t (observed) and t (model implied, bootstrap)"
  
  return(list(ks, c(tdist), c(as.logical(published))))
}