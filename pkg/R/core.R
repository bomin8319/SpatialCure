#' @useDynLib SpatialCure
#' @importFrom stats dgamma runif
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
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
betas.slice.sampling = function(Sigma.b, Y, X, betas, alpha, C, lambda, w, m, form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling(betas[p], p, Sigma.b, Y, X, betas, alpha, C, lambda, w, m, form = form)
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
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.betas.slice.sampling = function(betas.p, p, Sigma.b, Y, X, betas, alpha, C, lambda, w, m, lower = -Inf, upper = +Inf, form) {
  b0 = betas.p
  b.post0 = betas.post(b0, p, Sigma.b, Y, X, betas, alpha, C, lambda, form)
  
  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, X, betas, alpha, C, lambda, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, X, betas, alpha, C, lambda, form) <= b.post0) break
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
    b.post1 = betas.post(b1, p, Sigma.b, Y, X, betas, alpha, C, lambda, form)
    
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
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
gammas.slice.sampling = function(Sigma.g, Y, eXB, Z, gammas, C, lambda, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling(gammas[p], p, Sigma.g, Y, eXB, Z, gammas, C, lambda, w, m, form = form)
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
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.gammas.slice.sampling = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post(g0, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form)
  
  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form) <= g.post0) break
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
    g.post1 = gammas.post(g1, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form)
    
    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

#' @title lambda.slice.sampling
#' @description univariate slice sampling for lambda
#'
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#'
#' @return One sample update using slice sampling
#'
#' @export
lambda.slice.sampling = function(Y, eXB, alpha, C, lambda, w, m, lower = 0.01, upper = +Inf) {
  l0 = lambda
  l.post0 = lambda.post(Y, eXB, alpha, C, l0)
  
  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (lambda.post(Y, eXB, alpha, C, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (lambda.post(Y, eXB, alpha, C, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (lambda.post(Y, eXB, alpha, C, L) <= l.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (lambda.post(Y, eXB, alpha, C, R) <= l.post0) break
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
    l.post1 = lambda.post(Y, eXB, alpha, C, l1)
    
    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}


#' @title betas.post
#' @description log-posterior distribution of betas with pth element fixed as betas.p
#'
#' @param betas.p current value of the pth element of betas
#' @param p pth element
#' @param Sigma.b variance estimate of betas
#' @param Y response variable
#' @param X covariates for betas
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return log- posterior density of betas
#'
#' @export
betas.post = function(betas.p, p, Sigma.b, Y, X, betas, alpha, C, lambda, form) {
  betas[p] = betas.p
  #if (form %in% "Weibull") {
  #  eXB = exp(X %*% betas + 1 / lambda)
  #} else {
    eXB = exp(X %*% betas)
  #}
  lprior = dmvnorm(betas, rep(0, length(betas)), Sigma.b, log = TRUE)
  lpost = llikWeibull(Y, eXB, alpha, C, lambda) + lprior
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
#' @param lambda current value of lambda
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return log- posterior density of betas
#'
#' @export
gammas.post = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, lambda, form) {
  gammas[p] = gammas.p
  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas - 1/lambda))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  lprior = dmvnorm(gammas, rep(0, length(gammas)), Sigma.g, log = TRUE)
  lpost = llikWeibull(Y, eXB, alpha, C, lambda) + lprior
  return(lpost)
}

#' @title lambda.post
#' @description log-posterior distribution of lambda
#'
#' @param Y response variable
#' @param eXB exponentiated vector of covariates times betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param a shape parameter of gammas prior
#' @param b scale parameter of gammas prior
#'
#' @return log- posterior density of betas
#'
#' @export
lambda.post = function(Y, eXB, alpha, C, lambda, a = 1, b = 1) {
  lprior = dgamma(lambda, a, b, log = TRUE)
  lpost = llikWeibull(Y, eXB, alpha, C, lambda) + lprior
  return(lpost)
}

#' @title mcmcOF
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian parametric OF model
#'
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for betas
#' @param Z covariates for gammas
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param w size of the slice in the slice sampling for (betas, gammas, lambda)
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcOF <- function(Y, C, X, Z, N, burn, thin, w = c(1, 1, 1), m = 10, form) {
  p1 = dim(X)[2]
  p2 = dim(Z)[2]
  
  # initial values
  betas = rep(0, p1)
  gammas = rep(0, p2)
  lambda = 1
  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas - 1/lambda))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  Sigma.b = 10 * p1 * diag(p1)
  Sigma.g = 10 * p2 * diag(p2)
  betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gammas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  lambda.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    if (iter %% 1000 == 0) print(iter)
    if (iter > burn) {
    Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1 * diag(p1))
    Sigma.g = riwish(1 + p2, gammas %*% t(gammas) + p2 * diag(p2))
    }
    betas = betas.slice.sampling(Sigma.b, Y, X, betas, alpha, C, lambda, w[1], m, form = form)
    #if (form %in% "Weibull") {
    #  eXB = exp(X %*% betas + 1/lambda)
    #} else {
      eXB = exp(X %*% betas)
    #}
    gammas = gammas.slice.sampling(Sigma.g, Y, eXB, Z, gammas, C, lambda, w[2], m, form = form)
    if (form %in% "Weibull") {
      alpha = 1 / (1 + exp(-Z %*% gammas - 1/lambda))
    } else {
      alpha = 1 / (1 + exp(-Z %*% gammas))
    }
    if (form %in% "Weibull") {
      lambda = lambda.slice.sampling(Y, eXB, alpha, C, lambda, w[3], m)
    } 
    if (iter > burn & (iter - burn) %% thin == 0) {
      betas.samp[(iter - burn) / thin, ] = betas
      gammas.samp[(iter - burn) / thin, ] = gammas
      lambda.samp[(iter - burn) / thin] = lambda
    }
  }
  return(list(betas = betas.samp, gammas = gammas.samp, lambda = lambda.samp))
}

