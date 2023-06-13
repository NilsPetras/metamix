
# sig or nonsig only
LogLike <- Vectorize(
  function (d, t, nu, r, c, nonSigOnly){
    ncp <- d * r
    if(!nonSigOnly) {
      L  <-  sum(
        dt(t, nu, ncp,  log = TRUE) -
          pt(c, nu, ncp, lower.tail = FALSE,  log.p = TRUE))
    }

    if(nonSigOnly) {
      L  <-  sum(
        dt(t, nu, ncp,  log = TRUE) -
          pt(c, nu, ncp,  log.p = TRUE))
    }

    return(L)
  },  vectorize.args = "d")


# mixture

LogLikeMix <- Vectorize(
  function (d, l, t, nu, r, c, SigSuppress, TwoSided){
    p <- 1 / (1 + exp(-l))
    ncp <- d * r
    # for extreme values,  precision of dt()  /  pt() causes warnings
    if (!TwoSided) {
      MLE <- sum(log(MixPdf(t, p, nu, ncp, c, SigSuppress)))
    }
    if (TwoSided) {
      MLE <- sum(log(MixPdfTwoSided(t, p, nu, ncp, c)))
    }
    MLE
  },  vectorize.args = c("d",  "l"))

MixPdf <- Vectorize(
  function (t, p, nu, ncp, c, SigSuppress){
    beta <- pt(c, nu, ncp)
    if (!SigSuppress) {
      A <- 1 - p * beta
      I <- ifelse(t <= c,  1,  0)
    }

    if (SigSuppress) {
      A <- 1 - p * (1 - beta)
      ifelse(t <= c,  0,  1)
    }

    f <- dt(t, nu, ncp) * (1 - p * I) / A
    return(f)
  },  vectorize.args = c("t",  "nu",  "ncp", "c"))

MixPdfTwoSided<-Vectorize(
  function(t, p, nu, ncp, c){
    A <- dt(t, nu, ncp)
    B <- pt(c, nu, ncp)
    C <- pt(-c, nu, ncp)
    publish <- 1 - p * (B - C)
    I <- ifelse(abs(t) >= c, 0, 1)
    f <- A / publish * (1 - p * I)
    return(f)
  },  vectorize.args = c("t", "nu", "ncp", "c"))

LogLikeMixP0 <- Vectorize(
  function(d, t, nu, r){
    ncp <- d * r
    L <- sum(dt(t, nu, ncp,  log =  TRUE))
    L
  },  vectorize.args = c("d"))
