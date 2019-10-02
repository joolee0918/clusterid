#ifndef COMMONF_H
#define COMMONF_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

using namespace Rcpp;


IntegerVector findInterval(NumericVector x, NumericVector breaks);
double dAMH(NumericVector u, double rho, bool logf);
NumericVector vdAMH(NumericVector u1, NumericVector u2, NumericVector rho, bool logf);
double dClayton(NumericVector u, double rho, bool logf);
NumericVector vdClayton(NumericVector u1, NumericVector u2, double rho, bool logf);
double dFrank(NumericVector u, double rho, bool logf);
NumericVector vdFrank(NumericVector u1, NumericVector u2, double rho, bool logf);
double dcopf(NumericVector u,  double rho, bool logf, int copula);
NumericVector vdcopf(NumericVector u1, NumericVector u2,  double rho, bool logf, int copula);
double dccopf(NumericVector u, double u0, double rho, bool logf, int copula);
NumericVector vdccopf(NumericVector u1, NumericVector u2, NumericVector u0, double rho, bool logf, int copula);

double pClayton(NumericVector u, double rho);
double pAMH(NumericVector u, double rho);
double pccopf(NumericVector u, double u0, double rho, int copula);
double pcopf(NumericVector u, double rho, int copula);
double hf(NumericVector u, IntegerVector del, double rho, int copula);
double hcf(NumericVector u, IntegerVector del, double u0, double rho, int copula);
double hf_Clayton(NumericVector u, IntegerVector del, double rho);
double hf_Frank(NumericVector u, IntegerVector del, double theta);
double h1(double u1, double u2, double rho, int copula);
double hc1(double u1, double u2, double u0, double rho, int copula);
NumericVector vh1(NumericVector u1, NumericVector u2, double rho, int copula);
NumericVector vhc1(NumericVector u1, NumericVector u2, NumericVector u0, double rho, int copula);


NumericVector cumsum3(NumericVector x);
NumericVector cumprod3(NumericVector x);
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts,  int logf);
NumericVector Hpc(NumericVector x, NumericVector levels, NumericVector cuts,  int logf);
double ppc(double q, NumericVector levels, NumericVector cuts, int lower, int logf);
double dpc(double x, NumericVector levels, NumericVector cuts, int logf);
NumericVector vppc(NumericVector q, NumericVector levels, NumericVector cuts, int lower, int logf);
NumericVector vdpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf);
double pG0(arma::vec r_id, NumericVector G, double p);
double pG(arma::vec r_id, NumericVector G, double p);
double ff1(int j1, int j2, NumericVector vpx, NumericVector vpx2, double alpha,
           NumericVector lam01, double newrho, double rho, NumericVector exam_age, NumericVector cut_F,
           List LAM03, List LAM12, NumericVector cut1, NumericVector cut2, NumericVector IG,
           NumericVector SS, arma::vec rid, double pg0, double p, NumericVector w1,  NumericVector w2,
           NumericVector u1,  NumericVector u2,  NumericVector ww1,  NumericVector ww2,  NumericVector uu1,  NumericVector uu2);
#endif

