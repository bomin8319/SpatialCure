% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{W.post}
\alias{W.post}
\title{W.post}
\usage{
W.post(S, A, lambda, Y, X, W, betas, delta, C, rho)
}
\arguments{
\item{S}{spatial information (e.g. district)}

\item{A}{adjacency information corresponding to spatial information}

\item{lambda}{CAR parameter}

\item{Y}{response variable}

\item{X}{covariates for betas}

\item{W}{spatial random effects}

\item{betas}{current value of betas}

\item{delta}{probability of true censoring}

\item{C}{censoring indicator}

\item{rho}{current value of rho}
}
\value{
log- posterior density of betas
}
\description{
log-posterior distribution of W with sth element fixed as W.s
}
