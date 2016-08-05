#include <Rcpp.h>

// roxygen
//' Fisher information matrix for the two-parameter logistic (2PL) model.
//'
//' The mean of response variable is \deqn{f(x, \bold{\theta}) = \frac{1}{(1 + \exp(-b (x - a)))},}{f(x, \bold{\theta}) = 1/(1 + \exp(-b (x - a))),} where
//' \eqn{\bold{\theta} = (a, b)}.
//' @param x vector of design points. In IRT, \eqn{x} is the person ability parameter.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{\bold{\theta} = (a, b)}. In IRT parameter
//' \eqn{a} is the item difficulty parameter and parameter \eqn{b} is the item discrimination parameter.
//' @return Fisher information matrix.
//' @family FIM
//' @export
//' @details
//' There is no closed-form for the locally optimal design.
//'  For minimax and standardized D-optimal design, the optimal design is symmetric around point
//' \eqn{(a^L + a^U)/2}{(aL + aU)/2} where \eqn{a^L}{aL} and \eqn{a^U}{aU} are the
//' lower bound and upper bound for parameter \eqn{a}, respectively. In \code{\link{mica}},
//'  options \code{sym} and \code{sym_point} in \code{control} can be used to make the search
//'   for the optimal design easier.
//'
//' @importFrom Rcpp evalCpp
//' @useDynLib ICAOD
// [[Rcpp::export]]

Rcpp::NumericMatrix FIM_logistic(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
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
        constant = 1/(1+exp(-b*(x[i]-a)));
        constant = w[i]*constant*(1-constant);
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


