#include <RcppEigen.h>



// roxygen
//' Fisher information matrix for the four parameter logistic model.
//'
//' The mean of the response variable is
//'  \deqn{f(x, \bold{\theta}) = \frac{\theta_1}{1 + \exp(\theta_2 x + \theta_3)} + \theta_4,}{
//'  f(x, \bold{\theta})= \theta1/(1 + exp(\theta2*x + \theta3)) + \theta4,}
//'   where \eqn{\bold{\theta} = (\theta_1, \theta_2, \theta_3, \theta_4)}{\bold{\theta} = (\theta1, \theta2, \theta3, \theta4)}.
//' @param x vector of design points.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters
//'  \eqn{\bold{\theta} = (\theta_1, \theta_2, \theta_3, \theta_4)}{\bold{\theta} = (\theta1, \theta2, \theta3, \theta4)}.
//' @return Fisher information matrix.
//' @details The fisher information matrix does not depend on \eqn{\theta_4}{\theta4}.\cr
//' There is no analytical solution for the locally D-optimal design for this model.
//' @family FIM
//' @export
// [[Rcpp::export]]



Eigen::MatrixXd FIM_logistic_4par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  theta1, theta2,  theta3, theta4, constant1, constant2;
  theta1 = param[0];
  theta2 = param[1];
  theta3 = param[2];
  theta4 = param[3];

  theta4 = theta4 + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(4, 4);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < x.size(); i++)
  {
    constant1 = exp(x[i]*theta2 + theta3);
    constant2 = 1/(1 + constant1);


    Eigen::MatrixXd f(4, 1);
    f(0, 0) = constant2;
    f(1, 0) = pow(constant2, 2)*-theta1*x[i]*constant1;
    f(2, 0) = pow(constant2, 2)*-theta1*constant1;
    f(3, 0) = 1,


    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;

  }

  return Fisher_mat;
}
