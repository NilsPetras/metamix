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
#' @param ConCof confidence coefficient for calculating the confidence interval;
#'   .95 corresponds to a 95\% CI
#' @param TwoSided assumption: were the observed t-tests two-sided?
#' @param SigSuppress assumption: were significant results suppressed?
#'
#' @return list including the ks-test output, a matrix of the bootstrapped
#'   t-values, and a matrix of the bootstrapped decision to publish / not
#'   publish
tFit <- function(
    t, n1, n2, alpha, p, d, nrep, ConCof, TwoSided, SigSuppress
) {
  # still missing: if n2 = NULL one sample t test
  
  # non-mixture = mixture model with p_sp = 1 (already done in metamix)
  k <- length(n1)
  df <- n1 + n2 - 2
  
  # guess the number of needed repetitions to achieve nrep
  # (including d could improve this for d >> 0)
  # and repeat with a better guess if too low
  published <- matrix(1)
  while(min(colSums(published, na.rm = T)) < nrep) {
    if (exists("n_new")) {
      n <- n_new
    } else {
      n <- ceiling(nrep / (1 - p + (1 - ConCof)))
    }
    
    tdist <- matrix( # non-central, model implied distribution
      nrow = n,
      ncol = k)
    
    published <- matrix( # publication status
      stats::runif(n * k),
      nrow = n,
      ncol = k)
    
    counter <- matrix(
      nrow = n,
      ncol = k)
    
    for (j in 1:k) { # distribution differs by sample size
      tdist[ ,j] <- stats::rt(
        n,
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
      for (j in 1:k) {
        # always publish significant t-values
        published[tdist[,j] > tcrit[j],j] <- TRUE
        # if non-significant, selective reporting check decides
        published[tdist[,j] < tcrit[j],j] <- published[tdist[,j] < tcrit[j],j] > p
      }
    } else { # > and < swapped
      for (j in 1:k) {
        published[tdist[,j] < tcrit[j],j] <- TRUE
        published[tdist[,j] > tcrit[j],j] <- published[tdist[,j] > tcrit[j],j] > p
      }
    }
    published <- matrix(as.logical(published), ncol = k)
    
    n_new <- ceiling(1.01 * n * (nrep / min(colSums(published, na.rm = T))))
  }
  
  
  # trim to an equal number of published repetitions per sample size combination
  for (j in 1:k) {
    counter[published[ ,j],j] <- 1:sum(published[ ,j])
  }
  
  cutoff <- apply(counter, 2, function(x) {
    which(x == nrep)
  })
  
  for (j in 1:k) {
    tdistsave[(cutoff[j] + 1):n,j] <- NA
    published[(cutoff[j] + 1):n,j] <- NA
  }
  
  # kolmogoroff smirnov test of observed and bootstrapped distributions
  ks <- stats::ks.test(t, c(tdistsave[published]))
  
  # for print
  ks$data.name <- "t (observed) and t (model implied, bootstrap)"
  
  return(list(ks, tdistsave, published))
}