#include <Rcpp.h>

// roxygen
//' Fisher information matrix for the one-parameter logistic model (1PL or Rasch model).
//'
//' The mean of response variable is
//'  \deqn{f(x, a) = 1/(1 + \exp(-(x - a))).}{f(x, a) = \frac{1}{(1 + exp(-(x - a)))}.}
//'  This function returns Fisher information for the design \eqn{\xi} that is
//'  \deqn{M(\xi; a) = \sum_{i = 1}^kw_iM(x_i, a).}{M(\xi, a) = sum w_i M(x_i, a).}
//'   Here \eqn{M(x, a)}  is \eqn{g(x-a)}, where
//'  \eqn{g(z) = \frac{\exp(z)}{(1 + \exp(z))^2}}{g(z) = exp(z)/(1 + exp(z))^2}.
//'  denotes the standard logisitc density.
//'
//' @param x vector of design points. In IRT \code{x} is the person ability parameter.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param parameter \eqn{a}. In IRT, it is called difficulty parameter.
//' @return Fisher information as a one by one matrix.
//' @references
//' Grasshoff, U., Holling, H., & Schwabe, R. (2012). Optimal designs for the Rasch model. Psychometrika, 77(4), 710-723.
//' @details
//' The locally optimal design is a one point design with \eqn{x^* = a}{x* = a} and provides a value of
//' \eqn{M(\xi^*, a) =  1/4}{M(\xi*, a) =  1/4} for the information.
//' @family FIM
//' @export
// [[Rcpp::export]]


Rcpp::NumericMatrix FIM_logisitic_1par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    //Rcpp::Rcout<<"'x' and 'w' are not of the same length."<<std::endl;
    Rcpp::stop("'x' and 'w' are not of the same length.");

  }
  double a, z, constant = 0;
  a = param[0];

  size_t i;

  for(i=0; i < x.size(); i++)
  {
    z = x[i] - a;
    constant = exp(z)/pow(exp(z) + 1, 2);
    constant = w[i]*constant;
  }


  Rcpp::NumericMatrix Fisher(1, 1);
  Fisher(0, 0) = constant;

  return Fisher;
}


