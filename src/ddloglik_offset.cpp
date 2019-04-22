#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;

//' An Rcpp function that calculates the first and second order derivate and partial likelihood function.
//' @param N   The sample size
//' @param delta   Event indicator
//' @param z   Covariate matrix
//' @param beta   The coefficient
//' @param offset   The whole spanned matrix
//' 
//' @return the first and second derivate and partial likelihood
// [[Rcpp::export]]
List ddloglik(int n, arma::colvec delta, arma::mat z, arma::colvec beta, arma::mat offset){
  int col   = beta.n_rows;     // row number of beta

  arma::mat s0    = arma::flipud(arma::cumsum(arma::flipud(exp(z*beta+offset))));
  arma::mat s00   = arma::repmat(s0,1,col);      // span the matrix in block structure
  arma::mat s0_2  = arma::repmat(s0%s0,1,col);   // for calculating the hessian matrix
  arma::mat delta0= arma::repmat(delta,1,col);   // span delta to calculating the first derivative
  arma::mat bb    = arma::repmat(exp(z*beta+offset),1,col);  
  arma::mat s     = z%bb;                   
  arma::mat s1    = arma::flipud(arma::cumsum(arma::flipud(s)));
  arma::mat l1    = (z-s1/s00)%delta0;        
  arma::rowvec L1 = arma::sum(l1,0);            // get the first derivative
  arma::mat lambda= arma::cumsum(delta/s0);     // get the lambda
  arma::mat partial_likelihood = delta%((z*beta+offset)-log(s0)); // get the partial likelihood (according to definination of partial likelihood of cox model)
  arma::mat s2;         
  arma::mat l2;         
  arma::mat L2= arma::zeros<arma::mat>(col,col);// initilize Hessian matrix
  // Using a loop to calculate the Hessian Matrix
  for(int i = 0; i < col; i++)
  { s2      = s%(arma::repmat(z.col(i),1,col));
    s2      = arma::flipud(arma::cumsum(arma::flipud(s2)));
    arma::mat s11 = arma::repmat(s1.col(i),1,col);
    l2      = delta0%(s2/s00-s1%s11/(s0_2));
    L2.col(i) = sum(l2,0).t();
  }
  // return the result we need for calculation.
  List result;
  result["L1"]  = L1;
  result["L2"]  = L2;
  result["lambda"] = lambda;
  result["partial_likelihood"] = partial_likelihood;
  return result;
}








