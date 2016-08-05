//Dette and Biedermann (2003): Robust and Efficient designs for the Michaelis-Menten model.
// model = ax/(b+x) on [0, x0]
#include <Rcpp.h>

// roxygen
//' Fisher information matrix for the Michaelis-Menten model.
//'
//' The mean of the response variable is
//'  \deqn{f(x, \bold{\theta}) = \frac{ax}{(b + x)},}{f(x, \bold{\theta}) = ax/(b + x),}
//'  where \eqn{\bold{\theta} = (a, b)}.
//' @param x vector of design points.
//' @param w vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @param param vector of model parameters \eqn{\bold{\theta} = (a, b)}.
//' @return Fisher information matrix.
//' @references Rasch, D. (1990). Optimum experimental design in nonlinear regression. Communications in Statistics-Theory and Methods, 19(12), 4786-4806.
//' @details
//' There is an analytical solution for the locally D-optimal design. See Rasch (1990).
//' @family FIM
//' @export
// [[Rcpp::export]]


Rcpp::NumericMatrix FIM_michaelis(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
    if(x.size() != w.size()){
        //Rcpp::Rcout<<"'x' and 'w' are not of the same length."<<std::endl;
        Rcpp::stop("'x' and 'w' are not of the same length.");

    }
    double A = 0, B = 0, C = 0, a, b, constant, w_sum;
    a = param[0];
    b = param[1];

    size_t i;

    for(i=0; i < x.size(); i++)
    {
        constant =  pow(x[i], 2)/pow(b+x[i],2);
        //1/(1+exp(-b*(x[i]-a)));
        constant = w[i]*constant;
        A = constant + A;
        B = -a/(b+x[i]) * constant + B;
        C = pow(a, 2)/pow(b+x[i], 2) * constant + C;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = C;


    return Fisher_mat;
}


