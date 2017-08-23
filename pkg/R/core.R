#' @useDynLib SpatialCure
#' @importFrom stats dgamma runif rgamma dnorm
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack riwish
#' @importFrom coda mcmc
#' @importFrom mvtnorm rmvnorm dmvnorm
NULL


#' @title betas.slice.sampling
#' @description slice sampling for betas
#'
#' @param Sigma.b variance estimate of betas
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
betas.slice.sampling = function(Sigma.b, Y, X, W, betas, delta, C, rho, w, m, form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling(betas[p], p, Sigma.b, Y, X, W, betas, delta, C, rho, w, m, form = form)
  }
  return(betas)
}

#' @title univ.betas.slice.sampling
#' @description univariate slice sampling for betas.p
#'
#' @param betas.p current value of the pth element of betas
#' @param p pth element
#' @param Sigma.b variance estimate of betas
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.betas.slice.sampling = function(betas.p, p, Sigma.b, Y, X, W, betas, delta, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  b0 = betas.p
  b.post0 = betas.post(b0, p, Sigma.b, Y, X, W, betas, delta, C, rho, form)
  
  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, X, W, betas, delta, C, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, X, W, betas, delta, C, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, X, W, betas, delta, C, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, X, W, betas, delta, C, rho, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post(b1, p, Sigma.b, Y, X, W, betas, delta, C, rho, form)
    
    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}

#' @title gammas.slice.sampling
#' @description slice sampling for gammas
#'
#' @param Sigma.g variance estimate of gammas
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
gammas.slice.sampling = function(Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling(gammas[p], p, Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, form = form)
  }
  return(gammas)
}

#' @title univ.gammas.slice.sampling
#' @description univariate slice sampling for gammas.p
#'
#' @param gammas.p current value of the pth element of gammas
#' @param p pth element
#' @param Sigma.g variance estimate of gammas
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.gammas.slice.sampling = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post(g0, p, Sigma.g, Y, eXB, Z,gammas, C, rho, form)
  
  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post(g1, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form)
    
    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

#' @title rho.slice.sampling
#' @description univariate slice sampling for rho
#'
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#'
#' @return One sample update using slice sampling
#'
#' @export
rho.slice.sampling = function(Y, eXB, delta, C, rho, w, m, lower = 0.01, upper = +Inf) {
  l0 = rho
  l.post0 = rho.post(Y, eXB, delta, C, l0)
  
  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post(Y, eXB, delta, C, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (rho.post(Y, eXB, delta, C, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (rho.post(Y, eXB, delta, C, L) <= l.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (rho.post(Y, eXB, delta, C, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = rho.post(Y, eXB, delta, C, l1)
    
    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}

#' @title W.slice.sampling
#' @description slice sampling for W
#'
#' @param S spatial information (e.g. district)
#' @param A adjacency information corresponding to spatial information
#' @param lambda CAR parameter 
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
W.slice.sampling = function(S, A, lambda, Y, X, W, betas, delta, C, rho, w, m, form = form) {
  S_uniq = unique(cbind(S, W))
  W_uniq = S_uniq[order(S_uniq[,1]), 2]
  for (s in sample(1:nrow(A), nrow(A), replace = FALSE)) {
    W_uniq[s] = univ.W.slice.sampling(W_uniq[s], s, S, A, lambda, Y, X, W, betas, delta, C, rho, w, m, form = form)
  }
  W = W_uniq[S]
  return(W)
}

#' @title univ.W.slice.sampling
#' @description univariate slice sampling for W.s
#'
#' @param W.s current value of the sth element of W
#' @param s sth element
#' @param S spatial information (e.g. district)
#' @param A adjacency information corresponding to spatial information
#' @param lambda CAR parameter 
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.W.slice.sampling = function(W.s, s, S, A, lambda, Y, X, W, betas, delta, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  W0 = W.s
  W.post0 = W.post(W0, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form)
  
  u = runif(1, 0, w)
  L = W0 - u
  R = W0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (W.post(L, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form) <= W.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (W.post(R, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form) <= W.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (W.post(L, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form) <= W.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (W.post(R, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form) <= W.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat
  {
    W1 = runif(1, L, R)
    W.post1 = W.post(W1, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form)
    
    if (W.post1 >= W.post0) break
    if (W1 > W0) {
      R = W1
    } else {
      L = W1
    }
  }
  return(W1)
}

#' @title betas.post
#' @description log-posterior distribution of betas with pth element fixed as betas.p
#'
#' @param betas.p current value of the pth element of betas
#' @param p pth element
#' @param Sigma.b variance estimate of betas
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return log- posterior density of betas
#'
#' @export
betas.post = function(betas.p, p, Sigma.b, Y, X, W, betas, delta, C, rho, form) {
  betas[p] = betas.p
  eXB = exp(X %*% betas + W)
  lprior = dmvnorm(betas, rep(0, length(betas)), Sigma.b, log = TRUE)
  lpost = llikWeibull(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}

#' @title gammas.post
#' @description log-posterior distribution of gammas with pth element fixed as gammas.p
#'
#' @param gammas.p current value of the pth element of gammas
#' @param p pth element
#' @param Sigma.g variance estimate of gammas
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param rho current value of rho
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return log- posterior density of betas
#'
#' @export
gammas.post = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) {
  gammas[p] = gammas.p
  delta = 1 / (1 + exp(-Z %*% gammas))
  lprior = dmvnorm(gammas, rep(0, length(gammas)), Sigma.g, log = TRUE)
  lpost = llikWeibull(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}

#' @title W.post
#' @description log-posterior distribution of W with sth element fixed as W.s
#'
#' @param W.s current value of the sth element of W
#' @param s pth element
#' @param S spatial information (e.g. district)
#' @param A adjacency information corresponding to spatial information
#' @param lambda CAR parameter 
#' @param Y response variable
#' @param X covariates for betas
#' @param W spatial random effects
#' @param betas current value of betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return log- posterior density of betas
#'
#' @export
W.post = function(W.s, s, S, A, lambda, Y, X, W, betas, delta, C, rho, form) {
  W[which(S == s)] = W.s
  eXB = exp(X %*% betas + W)
  adj = which(A[s,] == 1)
  m_j = length(adj)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  W_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
  lprior = dnorm(W.s, W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  lpost = llikWeibull(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}


#' @title rho.post
#' @description log-posterior distribution of rho
#'
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param delta probability of true censoring
#' @param C censoring indicator
#' @param rho current value of rho
#' @param a shape parameter of gammas prior
#' @param b scale parameter of gammas prior
#'
#' @return log- posterior density of betas
#'
#' @export
rho.post = function(Y, eXB, delta, C, rho, a = 1, b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikWeibull(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}


#' @title lambda.gibbs.sampling
#' @description log-posterior distribution of rho
#'
#' @param S spatial information (e.g. district)
#' @param A adjacency information corresponding to spatial information
#' @param W spatial random effects
#' @param a shape parameter of gammas prior
#' @param b scale parameter of gammas prior
#'
#' @return log- posterior density of betas
#'
#' @export
lambda.gibbs.sampling <- function(S, A, W, a = 1, b = 1) {
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  J = nrow(S_uniq)
  sums = 0
  for (j in 1:J) {
    adj = which(A[j,]==1)
    m_j = length(adj)
    W_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
    W_j = S_uniq[j, 2]
    sums = sums + m_j/2 * (W_j-W_j_bar)^2
  }
  lambda = rgamma(1, J + a, sums + b)
  return(lambda)
}


#' @title mcmcSurv
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian survival (Weibull) model
#'
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for betas
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param w size of the slice in the slice sampling for (betas, gammas, rho)
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcSurv <- function(Y, C, X, N, burn, thin, w = c(1, 1, 1), m = 10, form) {
  p1 = dim(X)[2]
  # initial values
  betas = rep(0, p1)
  rho = 1
  W = rep(0, length(Y))
  delta = rep(1, length(Y))
  Sigma.b = 10 * p1 * diag(p1)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  rho.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    if (iter %% 500 == 0) print(iter)
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
    }
    betas = betas.slice.sampling(Sigma.b, Y, X, W, betas, delta, C, rho, w[1], m, form = form)
    eXB = exp(X %*% betas + W)
    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, eXB, delta, C, rho, w[3], m)
    } 
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      rho.samp[(iter - burn) / thin] = rho
    }
  }
  return(list(betas = betas.samp, rho = rho.samp))
}


#' @title mcmcCure
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian cure model
#'
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for betas
#' @param Z covariates for gammas
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param w size of the slice in the slice sampling for (betas, gammas, rho)
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcCure <- function(Y, C, X, Z, N, burn, thin, w = c(1, 1, 1), m = 10, form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  delta = 1 / (1 + exp(-Z %*% gammas))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    if (iter %% 500 == 0) print(iter)
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling(Sigma.b, Y, X, W, betas, delta, C, rho, w[1], m, form = form)
    eXB = exp(X %*% betas + W)
    gammas = gammas.slice.sampling(Sigma.g, Y, eXB, Z, gammas, C, rho, w[2], m, form = form)
    delta = 1 / (1 + exp(-Z %*% gammas))
    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, eXB, delta, C, rho, w[3], m)
    } 
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
    }
  }
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp))
}


#' @title mcmcSpatialCure
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian spatial cure model
#'
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for betas
#' @param Z covariates for gammas
#' @param S spatial information (e.g. district)
#' @param A adjacency information corresponding to spatial information
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param w size of the slice in the slice sampling for (betas, gammas, rho)
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcSpatialCure <- function(Y, C, X, Z, S, A, N, burn, thin, w = c(1, 1, 1), m = 10, form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  rho = 1
  lambda = 1
  W = rep(0, length(Y))
  delta = 1 / (1 + exp(- Z %*% gammas))
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  rho.samp = rep(NA, (N - burn) / thin)
  lambda.samp = rep(NA, (N - burn) / thin)
  W.samp = matrix(NA, nrow = (N - burn) / thin, ncol = nrow(A))
  for (iter in 1:N) {
    if (iter %% 500 == 0) print(iter)
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    #CAR model
    lambda = lambda.gibbs.sampling(S, A, W)
    W = W.slice.sampling(S, A, lambda, Y, X, W, betas, delta, C, rho, w[1], m, form = form)
    betas = betas.slice.sampling(Sigma.b, Y, X, W, betas, delta, C, rho, w[1], m, form = form)
    eXB = exp(X %*% betas + W)
    gammas = gammas.slice.sampling(Sigma.g, Y, eXB, Z, gammas, C, rho, w[2], m, form = form)
    delta = 1 / (1 + exp(- Z %*% gammas))
    
    if (form %in% "Weibull") {
      rho = rho.slice.sampling(Y, eXB, delta, C, rho, w[3], m)
    } 
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      rho.samp[(iter - burn) / thin] = rho
      lambda.samp[(iter - burn) / thin] = lambda
      S_uniq = unique(cbind(S, W))
      S_uniq = S_uniq[order(S_uniq[,1]),]
      W.samp[(iter - burn) / thin, ] = S_uniq[,2]
    }
  }
  return(list(betas = betas.samp, gammas = gammas.samp, rho = rho.samp, lambda = lambda.samp, W = W.samp))
}

