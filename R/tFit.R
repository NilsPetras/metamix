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
  # still missing: if n2 = NULL one sample t test
  
  # non-mixture = mixture model with p_sp = 1 (already done in metamix)
  
  df <- n1 + n2 - 2
  
  tdist <- matrix( # non-central, model implied distribution
    nrow = nrep,
    ncol = length(n1))
  
  published <- matrix( # publication status
    stats::runif(nrep * length(n1)),
    nrow = nrep,
    ncol = length(n1))
  
  for (j in seq_along(n1)) { # distribution differs by sample size
    tdist[ ,j] <- stats::rt(
      nrep,
      df[j],
      d * sqrt(n1[j] * n2[j] / (n1[j] + n2[j]))) # see Faul et al. (2007), G-Power paper, Table 3
  }
  tdistsave <- tdist
  
  if (TwoSided) {
    tcrit <- stats::qt(1 - alpha / 2, df)
    tdist <- abs(tdist)
  } else {
    tcrit <- stats::qt(1 - alpha, df)
  }
  
  if (!SigSuppress) {
      for (j in seq_along(n1)) {
        # always publish significant t-values
        published[tdist[,j] > tcrit[j],j] <- TRUE
        # if non-significant, selective reporting check decides
        published[tdist[,j] < tcrit[j],j] <- published[tdist[,j] < tcrit[j],j] > p
      }
  } else { # > and < swapped
      for (j in seq_along(n1)) {
        published[tdist[,j] < tcrit[j],j] <- TRUE
        published[tdist[,j] > tcrit[j],j] <- published[tdist[,j] > tcrit[j],j] > p
      }
  }
  
  # kolmogoroff smirnov test of observed and bootstrapped distributions
  ks <- stats::ks.test(t, c(tdistsave[as.logical(published)]))
  
  # for print
  ks$data.name <- "t (observed) and t (model implied, bootstrap)"
  
  return(list(ks, c(tdistsave), c(as.logical(published))))
}