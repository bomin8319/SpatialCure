% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{lambda.gibbs.sampling2}
\alias{lambda.gibbs.sampling2}
\title{lambda.gibbs.sampling2}
\usage{
lambda.gibbs.sampling2(S, A, W, V, a = 1, b = 1)
}
\arguments{
\item{S}{spatial information (e.g. district)}

\item{A}{adjacency information corresponding to spatial information}

\item{W}{spatial random effects}

\item{V}{spatial random effects}

\item{a}{shape parameter of gammas prior}

\item{b}{scale parameter of gammas prior}
}
\value{
log- posterior density of betas
}
\description{
log-posterior distribution of rho
}
