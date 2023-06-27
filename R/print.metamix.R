#' @export
print.metamix <- function(x, ...){
  
  est <- x$estimates
  estim <- data.frame(
    est = c(est$d_est),
    se = c(est$SE_d),
    ci.lower = c(est$CI_d[1]),
    ci.upper = c(est$CI_d[2]),
    chisq = c(est$X2),
    p.value = c(est$p))
  # in the mixture model there is a p_sp estimate
  if(!(x$model$SigOnly | x$model$nonSigOnly)){
    estim[2, ] <- c(est$p_est, est$SE_p, est$CI_p[1], est$CI_p[2], est$XX2, est$pp)
    row.names(estim)[2] <- "p_sp"
  }
  
  row.names(estim)[1] <- c("d")
  
  # actual printout
  
  # estimates: header
  cat("\n\t")
  if (!x$model$SigOnly & !x$model$nonSigOnly) cat("mixture")
  cat(" model of selective publishing\n\n")
  cat("assumed alpha =", x$model$alpha, "\n")
  if (x$model$TwoSided) cat("two-sided hypothesis tests\n")
  else cat("one-sided hypothesis tests\n")
  if (x$model$SigSuppress) cat("significant results suppressed\n\n")
  else cat("non-significant results suppressed\n\n")
  
  # estimates: estimates as data frame
  print(signif(estim, 3))
  
  # ks-test results (only changed the name of the variables)
  print(x$model_fit_test)
  
  invisible(x)
}
