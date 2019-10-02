
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>


using namespace Rcpp;

double log1mexp(double a)
{
  double result;
  if (a<log(2)) {
    result=log(-expm1(-a));
  }else{
    result=log1p(-exp(-a));
  }
  return result;
}



IntegerVector findInterval(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  for(it = x.begin(), out_it = out.begin(); it != x.end();
  ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return out;
}



double dAMH(NumericVector u, double rho, bool logf){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = 1;
  }else{
  res  = (-1 + pow(rho,2)*(-1 + u1 + u2 - u1*u2) - rho*(-2 + u1 + u2 + u1*u2)) / pow(-1 + rho*(-1 + u1)*(-1 + u2),3);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}


NumericVector vdAMH(NumericVector u1, NumericVector u2, NumericVector rho, bool logf){
  NumericVector res;
  int n = u1.size();
  int i;

  for(i=0; i< n; i++){
  if(rho[i] == 0) {
    res[i] = 1;
  } else{
    res[i]  = (-1 + pow(rho[i],2)*(-1 + u1[i] + u2[i] - u1[i]*u2[i]) - rho*(-2 + u1[i] + u2[i] + u1[i]*u2[i])) / pow(-1 + rho*(-1 + u1)*(-1 + u2),3);
  }
  if(logf == 1.0) res[i] = log(res[i]);
  }

  return(res);
}

double dClayton(NumericVector u, double rho, bool logf){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = 1;
  }else{
    res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}



NumericVector vdClayton(NumericVector u1, NumericVector u2, double rho, bool logf){
  NumericVector res;
  if(rho == 0) {
    res.fill(1);
  } else{
  res  = (1+rho)*pow(u1*u2, -rho-1)*pow((pow(u1, -rho) + pow(u2, -rho)-1), -1/rho-2);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}



double dFrank(NumericVector u, double rho, bool logf){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = 1;
  }else{
    res  = (rho*(exp(rho)-1.0)*exp(rho*u2+rho*u1+rho))/pow(exp(rho*u2+rho*u1)-exp(rho*u2+rho)-exp(rho*u1+rho)+exp(rho),2.0);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}


NumericVector vdFrank(NumericVector u1, NumericVector u2, double rho, bool logf){
  NumericVector res;
  if(rho == 0) {
    res.fill(1);
  } else{
    res  = (rho*(exp(rho)-1.0)*exp(rho*u2+rho*u1+rho))/pow(exp(rho*u2+rho*u1)-exp(rho*u2+rho)-exp(rho*u1+rho)+exp(rho),2.0);
  }
  if(logf == 1.0) res = log(res);
  return(res);
}

double dcopf(NumericVector u,  double rho, bool logf, int copula){
  double res=0;
  if(copula==1) {
    res = dClayton(u, rho, logf);
  }
  else res = dFrank(u, rho, logf);
}


NumericVector vdcopf(NumericVector u1, NumericVector u2,  double rho, bool logf, int copula){
  NumericVector res;
  if(copula==1) res = vdClayton(u1, u2, rho, logf);
  else res = vdFrank(u1, u2, rho, logf);
}

double dccopf(NumericVector u, double u0, double rho, bool logf, int copula){
  double res=0;
  if(copula==1) {
    res = dClayton(u, rho/(1+rho), logf);
  }
  else res = dAMH(u, 1-exp(-rho*u0), logf);
}

NumericVector vdccopf(NumericVector u1, NumericVector u2, NumericVector u0, double rho, bool logf, int copula){
  NumericVector res;
  if(copula==1) res = vdClayton(u1, u2, rho/(1+rho), logf);
  else res = vdAMH(u1, u2, 1-exp(-rho*u0), logf);
}

double pClayton(NumericVector u, double rho){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = u1*u2;
  } else{
  res = pow(pow(u1,-rho) + pow(u2,-rho)-1, -1/rho);
  }
  return(res);
}

double pFrank(NumericVector u, double rho){
  double t1, t2, res;

  if (rho>0) {
    t1=-log1p(exp(-rho) * expm1(rho-u[0]*rho)/expm1(-rho));
    t2=-log1p(exp(-rho) * expm1(rho-u[1]*rho)/expm1(-rho));
    res = -log1mexp(t1+t2-log1mexp(rho))/rho;
  } else {
    res =-1/rho * log(1 + exp(-(-log((exp(-rho * u[0]) - 1)/(exp(-rho) - 1)) + -log((exp(-rho * u[1]) - 1)/(exp(-rho) - 1)))) * (exp(-rho) - 1));
  }
}

double pAMH(NumericVector u, double rho){
  double res;
  double u1 = u[0];
  double u2 = u[1];
  if(rho==0) {
    res = u1*u2;
  } else{
    res =  u1 * u2 / (1 - rho * (1 - u1) * (1 - u2));
  }
  return(res);
}

double pccopf(NumericVector u, double u0, double rho, int copula){
  double res=0;
  if(copula==1) res = pClayton( u, rho/(1+rho));
  else res = pAMH(u, 1-exp(-u0*rho));
}


double pcopf(NumericVector u, double rho, int copula){
  double res=0;
  if(copula==1) res = pClayton( u, rho);
  else res = pFrank(u, rho);
}


double hf(NumericVector u, IntegerVector del, double rho, int copula){

  double res;
  int j;
  double term1, term2, term;

  if(copula==1){
    //Clayton
    term = 0;
    if(rho == 0){
      res =  pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);;
    }else{
      for(j=0; j<2; j++){
        term += pow(u[j], -rho);
      }

      term2 = 1;
      for(j=0; j<2; j++){
        term2 *= pow(u[j], del[j]);
      }

      res = pow(term2, -rho-1)*pow(term-2+1, -1/rho-1);
    }
  } else{
    //Frank
    if(rho == 0){
      res = pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);
    }else{

      res = -(exp(rho)*(exp(rho*u[0])*(1-del[0]) + exp(rho*u[1])*(1-del[1])-1.0))/(exp(rho*u[0]+rho*u[1])-exp(rho*u[1]+rho)-exp(rho*u[2]+rho)+exp(rho));
    }
  }

  return(res);
}

double hcf(NumericVector u, IntegerVector del, double u0, double rho, int copula){

  double res;
  int j;
  double term1, term2, term;
  double newrho;

  if(copula==1){
    //Clayton
    newrho = rho/(1+rho);
    term = 0;
    if(newrho == 0){
      res =  pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);;
    }else{
      for(j=0; j<2; j++){
        term += pow(u[j], -newrho);
      }

      term2 = 1;
      for(j=0; j<2; j++){
        term2 *= pow(u[j], del[j]);
      }

      res = pow(term2, -newrho-1)*pow(term-2+1, -1/newrho-1);
    }
  } else{
    //AMH
    newrho = 1-exp(-rho*u0);
    if(newrho == 0){
      res = pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);
    }else{
      res = (u[0]*(1-newrho*(1-u[0])))*(1-del[0])*(u[1]*(1-newrho*(1-u[1])))*(1-del[1])/pow(1-newrho*(1-u[0])*(1-u[1]), 2);
    }
  }

  return(res);
}

double hf_Clayton(NumericVector u, IntegerVector del, double rho){

  int j;

  double term1, term2, term;
  double res;

  term = 0;
  if(rho == 0){
    res =  pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);;
  }else{
  for(j=0; j<2; j++){
    term += pow(u[j], -rho);
  }

  term2 = 1;
  for(j=0; j<2; j++){
    term2 *= pow(u[j], del[j]);
  }

  res = pow(term2, -rho-1)*pow(term-2+1, -1/rho-1);
  }
  return(res);
}


double hf_Frank(NumericVector u, IntegerVector del, double theta){

  double res;

  if(theta == 0){
    res = pow(u[0], 1-del[0])*pow(u[1], 1-del[1]);
  }else{
    res = -(exp(theta)*(exp(theta*u[0])*(1-del[0]) + exp(theta*u[1])*(1-del[1])-1.0))/(exp(theta*u[0]+theta*u[1])-exp(theta*u[1]+theta)-exp(theta*u[2]+theta)+exp(theta));
  }
  return(res);
}


double h1(double u1, double u2, double rho, int copula){

  double res;

  if(rho == 0){
    res = u2;
  }else{
    if(copula == 1) {
     res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
    } else{
     res = -(exp(rho)*(exp(rho*u2)-1.0))/(exp(rho*u1+rho*u2)-exp(rho*u1+rho)-exp(rho*u2+rho)+exp(rho));
    }
  }
  return(res);
}


double hc1(double u1, double u2, double u0, double rho, int copula){

  double res;
  double newrho;
  if(rho == 0){
    res = u2;
  }else{
    if(copula == 1) {
      newrho = rho/(1+rho);
      res = pow(u1, -newrho-1)*pow((pow(u2, -newrho) + pow(u1, -newrho) -1), -1/newrho-1);
    } else{
      newrho = 1-exp(-rho*u0);
      res =  (u2*(1-newrho*(1-u2)))/pow(1-newrho*(1-u1)*(1-u2), 2);
    }
  }
  return(res);
}


NumericVector vh1(NumericVector u1, NumericVector u2, double rho, int copula){
  NumericVector res;
  if(rho == 0){
    res = u2;
  }else{
    if(copula == 1) {
      res = pow(u1, -rho-1)*pow((pow(u2, -rho) + pow(u1, -rho) -1), -1/rho-1);
    } else{
      res = -(exp(rho)*(exp(rho*u2)-1.0))/(exp(rho*u1+rho*u2)-exp(rho*u1+rho)-exp(rho*u2+rho)+exp(rho));
    }

  }
  return(res);
}

NumericVector vhc1(NumericVector u1, NumericVector u2, NumericVector u0, double rho, int copula){
  NumericVector res;

  if(rho == 0){
    res = u2;
  }else{
    if(copula == 1) {
      double newrho = rho/(1+rho);
      res = pow(u1, -newrho-1)*pow((pow(u2, -newrho) + pow(u1, -newrho) -1), -1/newrho-1);
    } else{
      NumericVector newrho = 1-exp(-rho*u0);
      res =  (u2*(1-newrho*(1-u2)))/pow(1-newrho*(1-u1)*(1-u2), 2);
    }

  }
  return(res);
}

NumericVector cumsum3(NumericVector x) {
  return cumsum(x); // compute + return result
}


NumericVector cumprod3(NumericVector x) {
  return cumprod(x); // compute + return result
}


//[[Rcpp::export()]]
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf)
{

 NumericVector cut = sort_unique(cuts);
 int p = levels.size();
 NumericVector y(x.size());
 cut.push_front(0);
 cut.push_back(R_PosInf);

  y[(cut[0] <= x) & (x < cut[1])] = levels[0];
      if (p > 1.5) {
        for (int i=1; i<p; i++) {
          y[(cut[i] <= x) & (x < cut[i + 1])] = levels[i];
        }
      }
      if (logf)
        y = log(y);
 return(y);
}

//[[Rcpp::export()]]
NumericVector Hpc(NumericVector x,  NumericVector levels, NumericVector cuts, int logf)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  LogicalVector who = (cut[0] <= x) & (x < cut[0 + 1]);
        if (sum(who)) {
          y[who] = x[who];
          y = y*levels[0];
        }
        double su = levels[0] * cut[1];
        if (p > 1.5) {
          for (int i = 1; i<p; i++) {
            who = (cut[i] <= x) & (x < cut[i + 1]);
            if (sum(who)) {
              NumericVector xwho= x[who];
              NumericVector tmpx = su + levels[i] * (xwho - cut[i]);
              y[who] = tmpx;
            }
            su = su + levels[i] * (cut[i + 1] - cut[i]);
          }
        }
        if (logf)
          y = log(y);
    return(y);
}


//[[Rcpp::export()]]
double ppc(double q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  double y;
  if (cuts[0]==0) {
     y = R::pexp(q, 1/levels[0], 0.0, 0.0);
    }else{
  NumericVector qq(1);
  qq[0] = q;
   y = Hpc(qq,  levels, cuts, 0.0)[0];
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
NumericVector vppc(NumericVector q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  NumericVector y(1);
  if (cuts[0]==0) {
    y = pexp(q, levels[0], lower, logf);
  }else{
    y = Hpc(q,  levels, cuts, 0.0);
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
double dpc(double x, NumericVector levels,  NumericVector cuts, int logf)
{

  double y;
  if (cuts[0]==0) {
    y = R::dexp(x, 1/levels[0], 0.0);
  }else{
    NumericVector xx(1);
    xx[0] = x;
  y = hpc(xx, levels, cuts,  0.0)[0] * ppc(x, levels, cuts, 0.0, 0.0);
  }
    if (logf)
      y = log(y);

    return(y);
}

//[[Rcpp::export()]]
NumericVector vdpc(NumericVector x, NumericVector levels,  NumericVector cuts, int logf)
{

  NumericVector y(x.size());
  if (cuts[0]==0) {
    y = dexp(x, levels[0]);
  }else{
    y = hpc(x, levels, cuts,  0.0) * vppc(x, levels, cuts, 0.0, 0.0);
  }
  if (logf)
    y = log(y);

  return(y);
}

//[[Rcpp::export()]]
NumericVector order_cpp(NumericVector x){
  match(x, clone(x).sort());
  return(x);
}

//[[Rcpp::export()]]
double pG0(arma::vec r_id, NumericVector G, double p){
  double q;
  double res = 0;

  q=1-p;
  //NumericVector order_id = match(r_id, clone(r_id).sort());
  arma::uvec order_rid = sort_index(r_id);
  NumericVector order_rid2 = as<NumericVector> (wrap(order_rid));

  G = G[order_rid2];

  if(sum(r_id)==2){
    if(G[0]==1 & G[1] ==1 ) res = pow(1-pow(q,2), 2);
    else if((G[0]==1 & G[1] ==0) | (G[0]==0 & G[1] ==1)) res = (1-pow(q,2))*pow(q,2);
    else res = pow(q,4);

  } else if(sum(r_id)==3){
    if(G[0]==1 & G[1] ==1 ) res = pow(p,2)*q+p;
      else if((G[0]==1 & G[1] ==0) | (G[0]==0 & G[1] ==1)) res = p*pow(q,2);
      else res = pow(q,3);

  } else if(sum(r_id)==4){
    if(G[0]==1 & G[1] ==1 ) res = 1.0/4.0*pow(p,2)*pow(1+p,2) + p*q*(2*p+1);
      else if((G[0]==1 & G[1] ==0) | (G[0]==0 & G[1] ==1)) res = 1.0/4.0*pow(p,2)*pow(q,2) + 1.0/2.0*p*pow(q,2)*(1+q);
        else res = 1.0/4.0*pow(q,2)*pow(1+q,2);
  }

  return(res);

}

//[[Rcpp::export()]]
double pG(arma::vec r_id, NumericVector G, double p){
  double q;
  double res = 0;
  q=1-p;

  arma::uvec order_rid = sort_index(r_id);
  NumericVector order_rid2 = as<NumericVector> (wrap(order_rid));
  G = G[order_rid2];
   if(sum(r_id) == 4){
    if(G[0]==1 & G[1] ==1 & G[2]==1) res = pow(p,2)*(1+2*q);
      else if(G[0]==1 & G[1]==1 & G[2]==0) res = pow(p,2)*pow(q,2);
      else if(G[0]==1 & G[1]==0 & G[2]==1) res = p*pow(q,2);
      else if(G[0]==1 & G[1]==0 & G[2]==0) res = p*pow(q,3);
      else if(G[0]==0 & G[1]==1 & G[2]==1) res = p*pow(q,2);
      else if(G[0]==0 & G[1]==1 & G[2]==0) res = p*pow(q,3);
      else if(G[0]==0 & G[1]==0 & G[2]==1) res = 0;
      else if(G[0]==0 & G[1]==0 & G[2]==0) res = pow(q,4);
  } else if(sum(r_id)==5){
    if(G[0]==1 & G[1] ==1 & G[2]==1) res = 1.0/4.0*pow(p,2)*(1+p)*(5-3*p) + 1/2*p*q*(p+p*q+1);
      else if(G[0]==1 & G[1]==1 & G[2]==0) res = 1.0/4.0*pow(p,2)*pow(q,2) + 1/2*p*pow(q,2);
      else if(G[0]==1 & G[1]==0 & G[2]==1) res = 1.0/4.0*pow(p,2)*pow(q,2) + 1/2*p*pow(q,2);
      else if(G[0]==1 & G[1]==0 & G[2]==0) res = 1.0/4.0*p*pow(q,2)*(1+q);
       else if(G[0]==0 & G[1]==1 & G[2]==1) res = 1.0/2.0*p*pow(q,2)*(1+p);
          else if(G[0]==0 & G[1]==1 & G[2]==0) res = 1.0/2.0*p*pow(q,3);
          else if(G[0]==0 & G[1]==0 & G[2]==1) res = 1.0/2.0*p*pow(q,3);
          else if(G[0]==0 & G[1]==0 & G[2]==0) res = 1.0/2.0*pow(q,3)*(1+q);
  } else if(sum(r_id)==6){
    if(G[0]==1 & G[1] ==1 & G[2]==1) res = 1.0/16.0*pow(p,2)*(1+3*p)*(7-3*p) + 1.0/4.0*p*q*(6*p + 3*p*q+2);
      else if(G[0]==1 & G[1]==1 & G[2]==0) res = 5.0/16.0*pow(p,2)*pow(q,2) + 1.0/4.0*p*pow(q,2)*(1+q);
        else if(G[0]==1 & G[1]==0 & G[2]==1) res = 5.0/16.0*pow(p,2)*pow(q,2) + 1.0/4.0*p*pow(q,2)*(1+q);
          else if(G[0]==1 & G[1]==0 & G[2]==0) res = 1.0/16.0*pow(p,2)*pow(q,2) + 1.0/8.0*p*pow(q,2)*(1+3*q);
            else if(G[0]==0 & G[1]==1 & G[2]==1) res = 5.0/16.0*pow(p,2)*pow(q,2) + 1.0/4.0*p*pow(q,2)*(1+q);
              else if(G[0]==0 & G[1]==1 & G[2]==0) res = 1.0/16.0*pow(p,2)*pow(q,2) + 1.0/8.0*p*pow(q,2)*(1+3*q);
                else if(G[0]==0 & G[1]==0 & G[2]==1) res = 1.0/16.0*pow(p,2)*pow(q,2) + 1.0/8.0*p*pow(q,2)*(1+3*q);
                  else if(G[0]==0 & G[1]==0 & G[2]==0) res = 1.0/16.0*pow(q,2)*pow(1+3*q,2);
  }
  return(res);
}



double ff1(int j1, int j2, double px, NumericVector vpx, NumericVector vpx2, double alpha,
           NumericVector lam01, double rho, NumericVector exam_age, NumericVector cut_F,
           List LAM03,  List LAM12, NumericVector cut1, NumericVector cut2, NumericVector IG,
           NumericVector SS, arma::vec rid, double pg0, double p, NumericVector w1,  NumericVector w2,
           NumericVector u1,  NumericVector u2,  NumericVector ww1,  NumericVector ww2,  NumericVector uu1,  NumericVector uu2, int copula) {

  double res;
  res = (pccopf(SS, px, rho, copula)*ppc(exam_age[j1], LAM03[j1], cut1, 0.0, 0.0)*ppc(exam_age[j2], LAM03[j2], cut2, 0.0, 0.0)
           + sum(w2*vhc1(vh1(vpx, vppc(u2, lam01*exp(alpha*IG[2]), cut_F, 0.0, 0.0), rho, copula), rep(SS[0], 20), vpx, rho, copula)
                   *vdcopf(vpx, vppc(u2, lam01*exp(alpha*IG[2]), cut_F, 0.0, 0.0), rho, 0.0, copula)*vdpc(u2, lam01*exp(alpha*IG[2]), cut_F, 0.0)*vppc(u2, LAM03[j2], cut2, 0.0, 0.0)
                   /vppc(u2, LAM12[j2], cut2, 0.0, 0.0))*ppc(exam_age[j2],LAM12[j2],  cut2, 0.0, 0.0)*ppc(exam_age[j1], LAM03[j1], cut1, 0.0, 0.0)
                   + sum(w1*vhc1(vh1(vpx, vppc(u1, lam01*exp(alpha*IG[1]), cut_F, 0.0, 0.0), rho, copula), rep(SS[1], 20), vpx, rho, copula)
                   *vdcopf(vpx, vppc(u1, lam01*exp(alpha*IG[1]), cut_F, 0.0, 0.0), rho, 0.0, copula)*vdpc(u1, lam01*exp(alpha*IG[1]), cut_F, 0.0)*vppc(u1, LAM03[j1], cut1, 0.0, 0.0)
                   /vppc(u1, LAM12[j1], cut1, 0.0, 0.0))*ppc(exam_age[j1],  LAM12[j1], cut1, 0.0, 0.0)*ppc(exam_age[j2], LAM03[j2], cut2, 0.0, 0.0)
                   +sum(ww1*ww2*vdccopf(vh1(vpx2, vppc(uu1, lam01*exp(alpha*IG[1]), cut_F, 0.0, 0.0), rho, copula) , vh1(vpx2, vppc(uu2, lam01*exp(alpha*IG[2]), cut_F, 0.0, 0.0), rho, copula),  vpx2, rho, 0.0, copula)
                   *vdcopf(vpx2, vppc(uu1, lam01*exp(alpha*IG[1]), cut_F, 0.0, 0.0), rho, 0.0, copula)*vdcopf(vpx2, vppc(uu2*exp(alpha*IG[2]), lam01, cut_F, 0.0, 0.0), rho, 0.0, copula)
                   *vdpc(uu1, lam01*exp(alpha*IG[1]), cut_F, 0.0)*vdpc(uu2, lam01*exp(alpha*IG[2]), cut_F, 0.0)*vppc(uu1, LAM03[j1], cut1, 0.0, 0.0)*ppc(exam_age[j1], LAM12[j1], cut1, 0.0, 0.0)
                   /vppc(uu1, LAM12[j1], cut1, 0.0, 0.0)*vppc(uu2, LAM03[j2],  cut2, 0.0, 0.0)
                   *ppc(exam_age[j2], LAM12[j2], cut2, 0.0, 0.0)/vppc(uu2, LAM12[j2], cut2, 0.0, 0.0)))*pG(rid, IG, p)/pg0;

  return(res);
}




