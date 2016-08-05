#include <Rcpp.h>

//#include <math.h>
//#include <vector>
//using namespace std;
//using namespace Rcpp;


//ON THE NUMBER OF SUPPORT POINTS OF MAXIMIN AND BAYESIAN OPTIMAL DESIGNS1
// BRAESS AND  DETTE (2007)
//a+exp(-bx)
//example3.3
// roxygen
//' Fisher information matrix for the two-parameter exponential model.
//'
//' The mean of response variable is \deqn{f(x, \bold{\theta}) = a + \exp(-bx)}.
//' @param x vector of design points.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{\bold{\theta} = (a, b)}.
//' @return Fisher information matrix.
//' @references Dette, H., & Neugebauer, H. M. (1997). Bayesian D-optimal designs for exponential regression models. Journal of Statistical Planning and Inference, 60(2), 331-349.
//' @family FIM
//' @details The Fisher information matrix does not depend on \code{a}.\cr
//' The locally D optimal design is independent of the nominal
//' value of \eqn{a} and is equally supported at \eqn{x = 0} and \eqn{x = 1/b}
//'  only when \eqn{x \in [0, 1]}{{x belongs to [0, 1]}}. See "Examples".
//' @export
//' @examples
//' \dontrun{
//' ### finding the locally optimal design for different values for design interval
//' mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = .0001, ux = 1, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 10, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = .0001, ux = 10, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' ## it seems for design interval x = [x_l, x_u], when x_l > 0,
//' ## the locally D-optimal design is a two-point equally weighted design
//' ## with x1 = x_l, x2 = x_u
//'
//' mica(fimfunc = "FIM_exp_2par", lx = .5, ux = 10, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = .0001, ux = 10, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = 1, ux = 10, lp = c(1, 2), up = c(1, 2),
//'         iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//'
//' mica(fimfunc = "FIM_exp_2par", lx = 2, ux = 10, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' mica(fimfunc = "FIM_exp_2par", lx = 3, ux = 9, lp = c(1, 2), up = c(1, 2),
//'      iter = 100, k = 2, type = "locally", control = list(seed = 215))
//'
//' }
// [[Rcpp::export]]




Rcpp::NumericMatrix FIM_exp_2par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }

  double A = 0, B = 0, C = 0, alpha, beta, constant, w_sum;
  alpha = param[0];
  beta = param[1];


  alpha = alpha + 0; //just to not get a warning


  size_t i;


    for(i=0; i < x.size(); i++)
    {
        constant = exp(-beta*x[i]);

        A = w[i] * 1 + A;
        B = w[i]*(-x[i]) * constant + B;
        C = w[i]*pow(x[i], 2)*exp(-2*beta*x[i]) + C;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = C;

    return Fisher_mat;
}
