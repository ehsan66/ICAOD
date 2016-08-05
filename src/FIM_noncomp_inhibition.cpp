#include <RcppEigen.h>





//Matlab code to produce the derivatives
//syms   I S V K_m K_ic
//noncomp_inhab = (V*S)/((K_m + S)*(1+I/K_ic));
//diff(noncomp_inhab, V)
//diff(noncomp_inhab, K_m)
//diff(noncomp_inhab, K_ic)


// roxygen
//' Fisher information matrix for the noncompetitive inhibition Michaelis-Menten model.
//'
//' The mean velocity of the reaction rate is
//'  \deqn{\eta = \frac{VS}{(K_m + S)(1 + \frac{I}{K_{ic}})}.}{\eta = VS/(Km + S)(1 + I/Kic).}
//'   Here, \eqn{S} is the substrate concentration, \eqn{I} is the inhibitor concentration,
//'    \eqn{V} is the maximum velocity of the enzyme,
//'     \eqn{K_{ic}}{Kic} is the dissociation constants and \eqn{K_m}{Km} is the Michaelis-Menten constant.
//'      Any design point is of the form \eqn{(S, I)}.
//' @param S vector of \code{S} component of design points. \code{S} is the substrate concentration.
//' @param I  vector of \code{I} component of design points. \code{I} is the inhibitor concentration.
//' @param w vector of corresponding weights for each design point. Its length must be equal to the length of \code{I} and \code{S}, and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{(V, K_m, K_{ic})}{(V, Km, Kic)}.
//' @return Fisher information matrix.
//' @references Bogacka, B., Patan, M., Johnson, P. J., Youdim, K., & Atkinson, A. C. (2011). Optimum design of experiments for enzyme inhibition kinetic models. Journal of biopharmaceutical statistics, 21(3), 555-572.
//' @details The model has an analytical solution for the locally D-optimal design. See Bogacka et al. (2011) for details.\cr
//' The optimal design does not depend on parameter \eqn{V}.
//' @family FIM
//' @export
// [[Rcpp::export]]


Eigen::MatrixXd FIM_noncomp_inhibition(const std::vector<double> S, const std::vector<double> I, const std::vector<double> w, const std::vector<double> param)
{
  if(I.size() != w.size() || S.size() != w.size()){
    Rcpp::stop("The length of 'I' or 'S' is not equal to the length of 'w'.");
  }
  double  V, K_m, K_ic, constant;
  V = param[0];
  K_m = param[1];
  K_ic = param[2];

  //a = a + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(3, 3);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < I.size(); i++)
  {

    Eigen::MatrixXd f(3, 1);
    constant = (K_m + S[i]);
    f(0, 0) = S[i]/((I[i]/K_ic + 1)*constant);
    f(1, 0) = -(S[i]*V)/((I[i]/K_ic + 1)*pow(constant,2));
    f(2, 0) = (I[i]*S[i]*V)/(pow(K_ic,2)*pow((I[i]/K_ic + 1),2)*constant);

    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;
  }

  return Fisher_mat;
}
