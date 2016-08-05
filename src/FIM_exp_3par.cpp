#include <RcppEigen.h>

//#include <math.h>
//#include <vector>
//using namespace std;
//using Eigen::MatrixXd;                  // variable size matrix, double  precision


// roxygen
//' Fisher information matrix for the three-parameter exponential model.
//'
//' The mean of response variable is
//'  \deqn{f(x, \bold{\theta}) = \theta_0 + \theta_1 \exp(\frac{x}{\theta_2})}{f(x, \bold{\theta}) = \theta0 + \theta1 \exp(x/\theta2)}.
//' @param x vector of design points.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters
//'  \eqn{\bold{\theta} = (\theta_0, \theta_1, \theta_2)}{\bold{\theta} =(\theta0, \theta1, \theta2)}.
//' @return Fisher information matrix.
//' @references Dette, H., Kiss, C., Bevanda, M., & Bretz, F. (2010). Optimal designs for the EMAX, log-linear and exponential models. Biometrika, 97(2), 513-518.
//' @family FIM
//' @details
//' The model has an analytical solution for the locally D-optimal design. See Dette et al. (2010) for more details.\cr
//' The Fisher information matrix does not depend on \eqn{\theta_0}{\theta0}.
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd FIM_exp_3par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  theta0, theta1,  theta2, constant;
  theta0 = param[0];
  theta1 = param[1];
  theta2 = param[2];

  theta0 = theta0 + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(3, 3);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < x.size(); i++)
  {

    constant = exp(x[i]/theta2);


    Eigen::MatrixXd f(3, 1);
    f(0, 0) = 1;
    f(1, 0) = constant;
    f(2, 0) = -theta1*x[i]*constant/pow(theta2, 2);



    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;

  }

  return Fisher_mat;
}
