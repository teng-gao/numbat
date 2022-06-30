// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppdbbinom
NumericVector cppdbbinom(const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& log_prob);
RcppExport SEXP _numbat_cppdbbinom(SEXP xSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP log_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_prob(log_probSEXP);
    rcpp_result_gen = Rcpp::wrap(cppdbbinom(x, size, alpha, beta, log_prob));
    return rcpp_result_gen;
END_RCPP
}
// logSumExp
double logSumExp(const arma::vec& x);
RcppExport SEXP _numbat_logSumExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp(x));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_compute
double likelihood_compute(Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m);
RcppExport SEXP _numbat_likelihood_compute(SEXP logphiSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type logphi(logphiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_compute(logphi, logprob, logPi, n, m));
    return rcpp_result_gen;
END_RCPP
}
// forward_backward_compute
Rcpp::NumericMatrix forward_backward_compute(Rcpp::NumericVector logphi, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m);
RcppExport SEXP _numbat_forward_backward_compute(SEXP logphiSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type logphi(logphiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_backward_compute(logphi, logprob, logPi, n, m));
    return rcpp_result_gen;
END_RCPP
}
// viterbi_compute
Rcpp::NumericVector viterbi_compute(Rcpp::NumericVector log_delta, Rcpp::NumericMatrix logprob, arma::cube logPi, int n, int m, Rcpp::NumericMatrix nu, Rcpp::NumericVector z);
RcppExport SEXP _numbat_viterbi_compute(SEXP log_deltaSEXP, SEXP logprobSEXP, SEXP logPiSEXP, SEXP nSEXP, SEXP mSEXP, SEXP nuSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_delta(log_deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type logprob(logprobSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi_compute(log_delta, logprob, logPi, n, m, nu, z));
    return rcpp_result_gen;
END_RCPP
}
// fit_lnpois_cpp
void fit_lnpois_cpp();
RcppExport SEXP _numbat_fit_lnpois_cpp() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    fit_lnpois_cpp();
    return R_NilValue;
END_RCPP
}
// poilog1
std::vector<double> poilog1(std::vector<int> x, std::vector<double> my, std::vector<double> sig);
RcppExport SEXP _numbat_poilog1(SEXP xSEXP, SEXP mySEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type my(mySEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(poilog1(x, my, sig));
    return rcpp_result_gen;
END_RCPP
}
// l_lnpois_cpp
std::vector<double> l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig);
RcppExport SEXP _numbat_l_lnpois_cpp(SEXP Y_obsSEXP, SEXP lambda_refSEXP, SEXP dSEXP, SEXP muSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type Y_obs(Y_obsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lambda_ref(lambda_refSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(l_lnpois_cpp(Y_obs, lambda_ref, d, mu, sig));
    return rcpp_result_gen;
END_RCPP
}
// allChildrenCPP
std::vector<std::vector<int>> allChildrenCPP(const arma::Mat<int> E);
RcppExport SEXP _numbat_allChildrenCPP(SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<int> >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(allChildrenCPP(E));
    return rcpp_result_gen;
END_RCPP
}
// CgetQ
arma::mat CgetQ(arma::mat logQ, std::vector<std::vector<int>> children_dict, arma::Col<int> node_order);
RcppExport SEXP _numbat_CgetQ(SEXP logQSEXP, SEXP children_dictSEXP, SEXP node_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type logQ(logQSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int>> >::type children_dict(children_dictSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type node_order(node_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(CgetQ(logQ, children_dict, node_order));
    return rcpp_result_gen;
END_RCPP
}
// score_tree_cpp
double score_tree_cpp(const arma::Mat<int> E, const arma::mat P);
RcppExport SEXP _numbat_score_tree_cpp(SEXP ESEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::Mat<int> >::type E(ESEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(score_tree_cpp(E, P));
    return rcpp_result_gen;
END_RCPP
}
// score_nni_parallel
NumericVector score_nni_parallel(List trees, arma::mat P);
RcppExport SEXP _numbat_score_nni_parallel(SEXP treesSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(score_nni_parallel(trees, P));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_numbat_cppdbbinom", (DL_FUNC) &_numbat_cppdbbinom, 5},
    {"_numbat_logSumExp", (DL_FUNC) &_numbat_logSumExp, 1},
    {"_numbat_likelihood_compute", (DL_FUNC) &_numbat_likelihood_compute, 5},
    {"_numbat_forward_backward_compute", (DL_FUNC) &_numbat_forward_backward_compute, 5},
    {"_numbat_viterbi_compute", (DL_FUNC) &_numbat_viterbi_compute, 7},
    {"_numbat_fit_lnpois_cpp", (DL_FUNC) &_numbat_fit_lnpois_cpp, 0},
    {"_numbat_poilog1", (DL_FUNC) &_numbat_poilog1, 3},
    {"_numbat_l_lnpois_cpp", (DL_FUNC) &_numbat_l_lnpois_cpp, 5},
    {"_numbat_allChildrenCPP", (DL_FUNC) &_numbat_allChildrenCPP, 1},
    {"_numbat_CgetQ", (DL_FUNC) &_numbat_CgetQ, 3},
    {"_numbat_score_tree_cpp", (DL_FUNC) &_numbat_score_tree_cpp, 2},
    {"_numbat_score_nni_parallel", (DL_FUNC) &_numbat_score_nni_parallel, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_numbat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
