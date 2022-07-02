#include <cmath>  // std::pow

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;

std::vector<double> l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig);

class fit_lnpois : public Functor {
public:
    double operator()(const arma::vec &x) override {

        std::vector<int> Y_obs = {10};
        std::vector<double> lambda_ref = {0.1};
        double d = 10;
        double mu = x[0];
        double sig = x[1];

        std::vector<double> l = l_lnpois_cpp(Y_obs, lambda_ref, d, mu, sig);

        return -log(l[0]);
    };
};

// [[Rcpp::export]]
void fit_lnpois_cpp() {
  fit_lnpois model;
  Roptim<fit_lnpois> opt("L-BFGS-B");
  opt.control.trace = 1;
  opt.set_hessian(true);
  arma::vec lower = {-arma::datum::inf, 0.01};
  opt.set_lower(lower);

  arma::vec x = {0, 1};
  opt.minimize(model, x);
  
  Rcpp::Rcout << "-------------------------" << std::endl;
  opt.print();
}