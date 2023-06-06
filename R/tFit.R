#' t Fit
#'
#' test the fit of the model against the observed distribution of t values
#' 
#' @param t
#' @param n1 
#' @param n2 
#' @param alpha 
#' @param p 
#' @param d 
#' @param nrep
#' @param TwoSided 
#' @param SigSuppress
#'
#' @return
#' @export
#'
#' @examples
tFit <- function(
    t, n1, n2, alpha, p, d, nrep, SigOnly, nonSigOnly, TwoSided, SigSuppress
) {
  # if n2 = NULL one sample t test
  
  tdist <- matrix(nrow = nrep, ncol = length(n1))
  
  df <- n1 + n2 - 2
  
    for (j in 1:length(n1)) {
      tdist[,j] <- rt(
        nrep,
        df[j],
        d * sqrt(n1[j] * n2[j] / (n1[j] + n2[j]))) # see Faul et al. (2007), G-Power paper, Table 3
    }
  
  published <- matrix(
    runif(nrep * length(n1)),
    nrow = nrep,
    ncol = length(n1))
  
  if (!SigSuppress) {
    if (!TwoSided) {
      tcrit <- qt(1 - alpha, df)
      for (j in 1:length(n1)) {
        published[tdist[,j] > tcrit[j],j] <- TRUE
        published[tdist[,j] < tcrit[j],j] <- published[tdist[,j] < tcrit[j],j] > p
      }
    } else {
      tcrit <- qt(1 - alpha / 2, df)
      for (j in 1:length(n1)) {
        published[abs(tdist[,j]) > tcrit[j],j] <- TRUE
        published[abs(tdist[,j]) < tcrit[j],j] <- published[abs(tdist[,j]) < tcrit[j],j] > p
      }
    }
    
  } else {
    if (!TwoSided) {
      tcrit <- qt(1 - alpha, df)
      for (j in 1:length(n1)) {
        published[tdist[,j] < tcrit[j],j] <- TRUE
        published[tdist[,j] > tcrit[j],j] <- published[tdist[,j] > tcrit[j],j] > p
      }
    } else {
      tcrit <- qt(1 - alpha / 2, df)
      for (j in 1:length(n1)) {
        published[abs(tdist[,j]) < tcrit[j],j] <- TRUE
        published[abs(tdist[,j]) > tcrit[j],j] <- published[abs(tdist[,j]) > tcrit[j],j] > p
      }
    }
    
  }

  
  ks <- ks.test(t, c(tdist[as.logical(published)]))
  
  return(list(ks, c(tdist), c(as.logical(published))))
}