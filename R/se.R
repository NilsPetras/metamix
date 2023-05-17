#' Standard error of d
#'
#' @param d
#' @param t
#' @param nu
#' @param r
#' @param c
#' @param ConCof
#' @param nonSigOnly
#'
#' @return
se <- function (d,t,nu,r,c,ConCof,nonSigOnly){
  h <- 0.0001   # differential
  g  <-  function(d) LogLike(d,t,nu,r,c,nonSigOnly)
  I  <-  (g(d+h) - 2*g(d) + g(d-h))/h^2 # Fisher info at x = d
  SE  <-  1/sqrt(-I) # Standard error of estimate
  p1 <- (1-ConCof)/2
  p2 <- ConCof+p1
  z <- qnorm(c(p1, p2),0,1)
  CI  <-  c(d+z[1]*SE, d+z[2]*SE)  # Confidence interval
  return(list(SE=SE, CI=CI))
}

#' standard error of d and p in mixed models
#'
#' @param d
#' @param p
#' @param t
#' @param nu
#' @param r
#' @param c
#'
#' @return
seMix <- function(d,p,t,nu,r,c,SigSuppress, TwoSided){
  h <- 0.0001
  if (p+h > 1)
    p <- p-h
  if (p-h < 0)
    p <- p+h

  I <- matrix(NA, 2,2)
  f <- function(d,p) LogLikeMix(d,log(p/(1-p)),t,nu,r,c,SigSuppress,TwoSided)
  I[1,1] <-  -(f(d+h,p)-2*f(d,p)+f(d-h,p))/h^2
  I[2,2] <-  -(f(d,p+h)-2*f(d,p)+f(d,p-h))/h^2
  I[1,2] <-  -(f(d+h,p+h)+f(d-h,p-h)-f(d-h,p+h)-f(d+h,p-h))/(4*h^2)
  I[2,1] <-  I[1,2]
  Cov <- solve(I)
  SE_d <- abs(sqrt(Cov[1,1]))
  SE_p <- abs(sqrt(Cov[2,2]))

  l <- log(p/(1-p))
  ff <- function(d,l) LogLikeMix(d,l,t,nu,r,c,SigSuppress,TwoSided)
  I[1,1] <-  -(ff(d+h,l)-2*ff(d,l)+ff(d-h,l))/h^2
  I[1,2] <-  -(ff(d+h,l+h)+ff(d-h,l-h)-ff(d-h,l+h)-ff(d+h,l-h))/(4*h^2)
  I[2,1] <-  I[1,2]
  I[2,2] <-  -(ff(d,l+h)-2*ff(d,l)+ff(d,l-h))/h^2
  Cov <- solve(I)
  SE_l <- abs(sqrt(Cov[2,2]))

  return(list(SE_d=SE_d, SE_p=SE_p, SE_l=SE_l))
}
