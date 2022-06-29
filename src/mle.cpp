#include <cmath>  // std::pow

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include <poilog.h>

using namespace roptim;

// [[Rcpp::export]]
std::vector<double> l_lnpois_cpp(std::vector<int> Y_obs, arma::vec lambda_ref, int d, double mu, double sig) {

    arma::vec mu_vec = mu + log(d * lambda_ref);

    return poilog1(Y_obs, arma::conv_to<std::vector<double>>::from(mu_vec), sig);
};

// class fit_lnpois : public Functor {
// public:
//     double operator()(const arma::vec &x) override {
//         double mu = x(0);
//         double sig = x(1);
//         arma::Col<int> Y_obs = {10};
//         arma::vec lambda_ref = {0.1};
//         double d = 10;

//         // return -l_lnpois_cpp(Y_obs, lambda_ref, d, mu, sig);
//     };
// };

// // [[Rcpp::export]]
// void fit_lnpois_cpp() {
//   fit_lnpois model;
//   Roptim<fit_lnpois> opt("BFGS");
//   opt.control.trace = 1;
//   opt.set_hessian(true);
  
//   arma::vec x = {0, 1};
//   opt.minimize(model, x);
  
//   Rcpp::Rcout << "-------------------------" << std::endl;
//   opt.print();
// }