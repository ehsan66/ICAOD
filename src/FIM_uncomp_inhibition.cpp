#include <RcppEigen.h>



//1.Optimum Design of Experiments for Enzyme Inhibition Kinetic Models
//Barbara Bogacka , Maciej Patan , Patrick J. Johnson , Kuresh Youdim , Anthony C. Atkinson
//Journal of Biopharmaceutical Statistics
// uncomp_inhab = V*S/(K_m + S*(1+I/K_iu));
// equation 5 and 7


//Matlab code to produce the derivatives

//syms   I S V K_m K_iu
//uncomp_inhab = V*S/(K_m + S*(1+I/K_iu));
//diff(uncomp_inhab, V)
//diff(uncomp_inhab, K_m)
//diff(uncomp_inhab, K_iu)


// roxygen
//' Fisher information matrix for the uncompetitive inhibition Michaelis-Menten model.
//'
//' The mean velocity of the reaction rate is
//'  \deqn{\eta =  \frac{VS}{(K_m + S(1 + \frac{I}{K_{iu}}))}.}{\eta =  VS/(Km + S(1 + I/Kiu)).}
//'  Here, \eqn{S} is the substrate concentration, \eqn{I} is the inhibitor concentration,
//'   \eqn{V} is the maximum velocity of the enzyme, \eqn{K_{iu}}{Kiu} is the dissociation constants and
//'    \eqn{K_m}{Km} is the Michaelis-Menten constant. Any design point is of the form \eqn{(S, I)}.
//' @param S vector of \code{S} component of design points. \code{S} is the substrate concentration.
//' @param I  vector of \code{I} component of design points. \code{I} is the inhibitor concentration.
//' @param w vector of corresponding weights for each design point.
//'  Its length must be equal to the length of \code{I} and \code{S}, and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{(V, K_m, K_{iu})}{(V, Km, Kiu)}.
//' @return Fisher information matrix.
//' @references Bogacka, B., Patan, M., Johnson, P. J., Youdim, K., & Atkinson, A. C. (2011). Optimum design of experiments for enzyme inhibition kinetic models. Journal of biopharmaceutical statistics, 21(3), 555-572.
//' @details The model has an analytical solution for the locally D-optimal design. See Bogacka et al. (2011) for details.\cr
//' The optimal design does not depend on parameter \eqn{V}.
//' @family FIM
//' @export
// [[Rcpp::export]]


Eigen::MatrixXd FIM_uncomp_inhibition(const std::vector<double> S, const std::vector<double> I, const std::vector<double> w, const std::vector<double> param)
{
  if(I.size() != w.size() || S.size() != w.size()){
    Rcpp::stop("The length of 'I' or 'S' is not equal to the length of 'w'.");
  }
  double  V, K_m, K_iu, constant;
  V = param[0];
  K_m = param[1];
  K_iu = param[2];

  //a = a + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(3, 3);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < I.size(); i++)
  {

    Eigen::MatrixXd f(3, 1);
    constant = (I[i]/K_iu + 1);
    f(0, 0) = S[i]/(K_m + S[i]*constant);
    f(1, 0) = -(S[i]*V)/pow((K_m + S[i]*constant),2);
    f(2, 0) = (I[i]*pow(S[i],2)*V)/(pow(K_iu,2)*pow((K_m + S[i]*constant),2));

    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;
  }

  return Fisher_mat;
}
