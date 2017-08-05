#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp; 


// **********************************************************//
//     	           Likelihood function            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull (arma::vec Y,
                    arma::vec eXB, 
                    arma::vec alpha,
                    arma::vec C,
                    double lambda) {
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
  arma::vec dexp2 = pow(eXB % Y, lambda - 1);
  arma::vec dexp3 = pow(eXB % Y, lambda);
  arma::vec llik1 = (1 - alpha) % dexp1 + lambda * alpha % eXB % dexp2 % dexp1;
  arma::vec lalpha = log(1 - alpha);
  
  arma::uvec ids0 = find(llik1 == 0);
  llik1.elem(ids0).fill(exp(-740));
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp3 == arma::datum::inf);
  dexp3.elem(ids2).fill(exp(700));
  arma::uvec ids3 = find(lalpha == -arma::datum::inf);
  lalpha.elem(ids3).fill(-740);
  arma::vec llik = C % log(llik1) + (1 - C) % (lalpha - dexp3);
  return sum(llik);
}


// **********************************************************//
//     	             Likelihood function                     //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull2(arma::vec Y,
					          arma::vec eXB, 
				          	arma::vec alpha,
				          	arma::vec C,
					          double lambda) {
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
	arma::vec dexp2 = pow(eXB % Y, lambda - 1);
	arma::vec dexp3 = pow(eXB % Y, lambda);
	arma::vec llik1 = (1 - alpha) + lambda * alpha % eXB % dexp2 % dexp1;
	arma::uvec ids1 = find(llik1 == arma::datum::inf);
	llik1.elem(ids1).fill(exp(700));
	arma::uvec ids0 = find(llik1 == 0);
	llik1.elem(ids0).fill(exp(-740));
	arma::uvec ids2 = find(dexp3 == arma::datum::inf);
	dexp3.elem(ids2).fill(exp(700));
	arma::vec llik = C % log(llik1) + (1 - C) % (log(alpha) - dexp3);
	return sum(llik);
}

