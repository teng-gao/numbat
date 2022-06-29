#include <RcppArmadillo.h>
#pragma once
// Function declarations
double poilog(int x, double my, double sig);
double maxf(int x, double my, double sig);
double upper(int x, double m, double my, double sig);
double lower(int x, double m, double my, double sig);
double my_f(double z, int x, double my, double sig, double fac);
void my_f_vec(double *z, int n, void *p);
std::vector<double> poilog1(std::vector<int> x, arma::vec my, double sig);
