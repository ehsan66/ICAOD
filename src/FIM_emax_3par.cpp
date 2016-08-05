#include <RcppEigen.h>



// roxygen
//' Fisher information matrix for the three-parameter emax model.
//'
//' The mean of response variable is
//'  \deqn{f(x, \bold{\theta}) = \theta_0 + \frac{\theta_1 x}{(x + \theta_2)}}{f(x, \bold{\theta}) = \theta0 + \theta1 x/(x + \theta2)}.
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


Eigen::MatrixXd FIM_emax_3par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
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

    constant = x[i] + theta2;


    Eigen::MatrixXd f(3, 1);
    f(0, 0) = 1;
    f(1, 0) = x[i]/constant;
    f(2, 0) = -theta1*x[i]/pow(constant, 2);



    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;

  }

  return Fisher_mat;
}
