#include <Rcpp.h>


// roxygen
//' Fisher information matrix for the power logistic model.
//'
//' The mean of response variable is
//'  \deqn{f(x; a, b, s) =  \frac{1}{(1 + \exp(-b (x - a)))^s},}{f(x; a, b, s) =  1/(1 + \exp(-b (x - a)))^s,}
//' @param x vector of design points.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{(a, b)}.
//' @param s power parameter.
//' @return Fisher information matrix.
//' @details
//' There is no analytical solution for the locally D-optimal design. Parameter \eqn{s} must be
//' passed by \code{...} in most of the functions like \code{\link{mica}}.
//' @family FIM
//' @export
// [[Rcpp::export]]


Rcpp::NumericMatrix FIM_power_logistic(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param, const double s)
{
    if(x.size() != w.size()){
      Rcpp::stop("'x' and 'w' are not of the same length.");

    }
    double A = 0, B = 0, C = 0, a, b, constant, w_sum;
    a = param[0];
    b = param[1];

    size_t i;

    for(i=0; i < x.size(); i++)
    {
        constant = 1/pow((1+exp(-b*(x[i]-a))), s);
        constant = w[i]*pow(s, 2)*constant*pow(1-pow(constant, 1/s),2)/(1-constant);
        A = pow(b,2)* constant + A;
        B = -b * (x[i] - a) * constant + B;
        C = pow((x[i] - a), 2) * constant + C;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = C;

    return Fisher_mat;
}


