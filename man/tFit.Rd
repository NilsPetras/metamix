% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tFit.R
\name{tFit}
\alias{tFit}
\title{t Fit}
\usage{
tFit(t, n1, n2, alpha, p, d, nrep, ConCof, TwoSided, SigSuppress)
}
\arguments{
\item{t}{observed t-values}

\item{n1}{observed sample size of group 1}

\item{n2}{observed sample size of group 2}

\item{alpha}{assumed alpha of observed t-tests}

\item{p}{ml estimate of the probability of selective publishing}

\item{d}{ml estimate of the true effect size}

\item{nrep}{number of bootstrap repetitions (total number will be nrep *
length(n1))}

\item{ConCof}{confidence coefficient for calculating the confidence interval;
.95 corresponds to a 95\% CI}

\item{TwoSided}{assumption: were the observed t-tests two-sided?}

\item{SigSuppress}{assumption: were significant results suppressed?}
}
\value{
list including the ks-test output, a matrix of the bootstrapped
  t-values, and a matrix of the bootstrapped decision to publish / not
  publish
}
\description{
Test the fit of the model against the observed distribution of t-values.
}
