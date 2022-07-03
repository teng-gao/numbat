#include <cmath>  // std::pow

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace roptim;

#include <RcppParallel.h>
using namespace RcppParallel;

double l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig);

class fit_lnpois : public Functor {

    public:
    
    const std::vector<int> Y_obs;
    
    const std::vector<double> lambda_ref;

    const int d;

    // initialize with source and destination
    fit_lnpois(std::vector<int> Y_obs, const std::vector<double> lambda_ref, const int d): 
        Y_obs(Y_obs), lambda_ref(lambda_ref), d(d) {}

    double operator()(const arma::vec &x) override {
        return -l_lnpois_cpp(Y_obs, lambda_ref, d, x[0], x[1]);
    };
};

// [[Rcpp::export]]
arma::rowvec fit_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d) {

  fit_lnpois model(Y_obs, lambda_ref, d);

  Roptim<fit_lnpois> opt("L-BFGS-B");
  opt.control.trace = 0;
  opt.set_hessian(false);
  arma::vec lower = {-arma::datum::inf, 0.01};
  opt.set_lower(lower);

  arma::vec x = {0, 1};
  opt.minimize(model, x);

  return opt.par().t();
}

// struct fit_worker : public Worker {

//     const arma::Mat<int> count_mat;
//     const std::vector<double> lambda_ref;
//     RMatrix<double> params;

//     fit_worker(const arma::Mat<int> count_mat, std::vector<double> lambda_ref, Rcpp::NumericMatrix params): 
//         count_mat(count_mat), lambda_ref(lambda_ref), params(params) {}

//     void operator()(std::size_t begin, std::size_t end) {
//         for (std::size_t i = begin; i < end; i++) {
//             arma::Col<int> counts = count_mat.col(i);
//             int d = sum(counts);
//             std::vector<int> counts_vec = arma::conv_to< std::vector<int> >::from(counts);
//             arma::rowvec res = fit_lnpois_cpp(counts_vec, lambda_ref, d);
//             params(i,0) = res(0);
//             params(i,1) = res(1);
//         }
//     }
// };

// // [[Rcpp::export]]
// Rcpp::NumericMatrix fit_lnpois_parallel(arma::Mat<int> count_mat, std::vector<double> lambda_ref) {
    
//     int n = count_mat.n_cols;

//     Rcpp::NumericMatrix params(n,2);

//     fit_worker fit_worker(count_mat, lambda_ref, params);

//     parallelFor(0, n, fit_worker);

//     return params;

// }