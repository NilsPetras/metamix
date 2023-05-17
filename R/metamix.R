#' Meta-analysis Mixture Model
#'
#' Estimate the true meta-analytical effect size corrected for publication bias.
#'
#' @param t vector of t-values
#' @param n1 vector of sample sizes
#' @param n2 second vector of sample sizes in case of two-sample t-tests
#' @param SigOnly did only significant study results enter the meta-analysis?
#' @param nonSigOnly did only non-significant study results enter the
#'   meta-analysis?
#' @param alpha significance level (if one-sided: right tale)
#' @param ConCof confidence coefficient for calculating the confidence interval;
#'   .95 corresponds to a 95% CI
#' @param TwoSided is the t-test two-sided?
#' @param SigSuppress if both significant and non-significant study results
#'   entered the meta-analysis, were significant results suppressed?
#'
#' @return
#' @export
#'
#' @examples
metamix <- function(
    t,
    n1, n2 = NULL,
    SigOnly = FALSE,
    nonSigOnly = FALSE,
    alpha = .05,
    ConCof = .95,
    TwoSided = FALSE,
    SigSuppress = FALSE) {

  if (SigOnly | nonSigOnly) {
    y <- EstimationSigOnly(
      t = t,
      n1 = n1,
      n2 = n2,
      alpha = alpha,
      ConCof = ConCof,
      nonSigOnly = nonSigOnly)
  } else {
    y <- EstimationMix(
      t = t,
      n1 = n1,
      n2 = n2,
      alpha = alpha,
      ConCof = ConCof,
      SigSuppress = SigSuppress,
      TwoSided = TwoSided
    )
  }

  return(y)
}

EstimationSigOnly <- function(
    t,
    n1, n2,
    alpha,
    ConCof,
    nonSigOnly){

  #### Output
  #  d_est: estimate of delta
  #  SE: standard error of estimate
  #  CI: confidence interval of delta
  #  X2: square of likelihood ratio test, df=1, H0: d=0
  #  p: p-value associated with X2

  # Starting values of d for nlminb-optimization
  dset <- c(0.6, 0.3) # removed 0 as the third starting value because it hung up estimation for example (Table 1 in Ulrich et al., 2018)

  # Compute critical values for all studies

  if(!is.null(n2)){
    c <-  qt(1-alpha, n1+n2-2)
    nu <- n1+n2-2
    r <- sqrt(n1*n2/(n1+n2))
  }

  if(is.null(n2)){
    c <-  qt(1-alpha, n1-1)
    nu <- n1-1
    r <- sqrt(n1)
  }

  f  <-  function(d)
    -LogLike(d,t,nu,r,c,nonSigOnly)

  # Convergence check with different starting values
  j <- 0
  fval <- dmax <- rep(NA, length(dset))
  for(d0 in dset){
    j <- j+1
    tmp  <-  nlminb(d0, f)
    fval[j] <- tmp$objective
    dmax[j] <- tmp$par
  }

  y <- min(fval, na.rm = TRUE)
  i <- which.min(fval)

  if (abs(var(dmax)) > 1e-3){
    cat('Check starting value \n')
    print(list(t=t, t.mean=mean(t),
               d_start=dset,
               dmax=dmax))
    d_est=NA; SE=NA
    CI <- c(NA, NA)
    X2=NA; p=NA
  }

  d_est  <-  dmax[i]
  SE_CI  <-  se(d_est,t,nu,r,c,ConCof,nonSigOnly)
  X2  <-  abs(-2*(LogLike(0,t,nu,r,c,nonSigOnly)-LogLike(d_est,t,nu,r,c,nonSigOnly)))

  return(list( d_est=d_est,
        SE=SE_CI$SE,
        CI=SE_CI$CI,
        X2=X2,
        p=1-pchisq(X2,1)))
}

EstimationMix <- function(
    t,
    n1,n2,
    alpha,
    ConCof,
    SigSuppress,
    TwoSided){

  # mixture model, two-sample t-test

  #### Output
  # d_est:  estimate of delta
  # SE_d:   standard error of d_est
  # CI_d:   confidence interval for delta
  # p_est:  estimate of p
  # SE_est: standard error of p_est
  # CI_p:   confidence interval for p
  # X2:     Chi-square of likelihood ratio test, df=1, H0: d=0
  # p:      p-value associated with X2
  # XX2:    Chi-square of likelihood ratio test, df=1, H0: p=0
  # pp:     p-value associated with XX2

  ##  Define starting values of p and d for nlminb
  pstart <-  c(.2, .5, .8)
  dstart <- c(0, .2, .5)

  ## Compute critical values for all studies;
  #  r is needed to compute delta
  if (!is.null(n2)) {
    nu <- n1+n2-2
    c <- qt(1-alpha, nu)
    r <- sqrt(n1*n2/(n1+n2))
  }

  if(is.null(n2)) {
    nu <- n1-1
    c <- qt(1-alpha, nu)
    r <- sqrt(n1)
  }


  ##  Estimate d and p
  f1  <-  function(x)
    -LogLikeMix(x[1],x[2],t,nu,r,c,SigSuppress,TwoSided)

  # Checks local minima for different starting values defined above
  fmax <- Inf
  i <- 0
  fs <- rep(NA, length(dstart))
  for (p0 in pstart){
    for (d0 in dstart){
      i <- i+1
      # many warnings due to precision (does not affect results)
      suppressWarnings(
        tmp  <-  nlminb(c(d0, log(p0/(1-p0))), f1)
      )
      fs[i] <- tmp$objective

      if (tmp$objective<fmax){
        d_est <- tmp$par[1]
        l_est <- tmp$par[2]
        fmax <- tmp$objective
      }
    }
  }
  if (sd(fs)>1e-3){
    cat('Check starting values\n')
    print(fs)
  }

  p_est <- 1/(1+exp(-l_est))
  SEtmp <- seMix(d_est,p_est,t,nu,r,c,SigSuppress, TwoSided)
  SE_d <- SEtmp$SE_d
  SE_l <- SEtmp$SE_l
  SE_p <- SEtmp$SE_p

  ## Confidence Intervals
  p1 <- (1-ConCof)/2
  p2 <- ConCof+p1
  z <- qnorm(c(p1, p2),0,1)
  CI_d  <-  d_est+z*SE_d
  CI_l  <-  l_est+z*SE_l
  CI_p  <-  1/(1+exp(-CI_l))

  ## Compute likelihood ratio test for d
  f0 <- function(y) -LogLikeMix(0,y,t,nu,r,c,SigSuppress,TwoSided)
  y <- nlminb(log(p_est/(1-p_est)),f0)$par
  p_0 <- 1/(1+exp(-y))
  X2 <- -2*(LogLikeMix(0,log(p_0/(1-p_0)),t,nu,r,c,SigSuppress,TwoSided) -
              LogLikeMix(d_est,log(p_est/(1-p_est)),t,nu,r,c,SigSuppress,TwoSided))
  p <- 1-pchisq(abs(X2),1)

  ## Compute likelihood ratio test for p
  f0 <- function(y) -LogLikeMixP0(y,t,nu,r)
  y <- nlminb(d_est,f0)$par
  XX2 <- -2*(LogLikeMixP0(y,t,nu,r) -
               LogLikeMix(d_est,log(p_est/(1-p_est)),t,nu,r,c,SigSuppress,TwoSided))
  pp <- 1-pchisq(abs(XX2),1)

  list(d_est=d_est, SE_d=SE_d, CI_d=CI_d,
       p_est=p_est, SE_p=SE_p, CI_p=CI_p,
       X2=X2, p=p, XX2=XX2, pp=pp)
}

