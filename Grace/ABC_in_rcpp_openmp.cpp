#include <Rcpp.h>
#include <RcppParallel.h>
#include <omp.h>
#include <sitmo.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(sitmo)]]
// [[Rcpp::plugins(openmp)]]

// Function to simulate data from SIR model
// [[Rcpp::export]]
NumericVector SIR(double s0, double i0, double r0, double beta, double gamma) {
  NumericVector s(20);
  NumericVector i(20);
  NumericVector r(20);
  
  s[0] = s0;
  i[0] = i0;
  r[0] = r0;

  for (int t = 1; t < s.size(); t++) {
    s[t] = s[t-1] - beta * s[t-1] * i[t-1];
    i[t] = i[t-1] + beta * s[t-1] * i[t-1] - gamma * i[t-1];
    r[t] = r[t-1] + gamma * i[t-1];
  }

  return i;
}

// [[Rcpp::export]]
double unif_sitmo(int seed) {
  uint32_t coreseed = static_cast<uint32_t>(seed);
  sitmo::prng eng(coreseed);
  double mx = sitmo::prng::max();
  double x = eng() / mx;
  return x;
}

double calc_dist(NumericVector x_sim, NumericVector x) {
  double total = 0.0;
  #pragma omp parallel reduction( + : total)
  {
    double inner_sum = 0.0;
    #pragma omp for
    for (int i=0; i<x_sim.size(); i++) {
      inner_sum += pow(x[i]-x_sim[i], 2);
    }
    total += inner_sum;
  }
  return total;
}

// [[Rcpp::export]]
NumericMatrix ABC(int n, double eps, int p, NumericVector x, int ncores)
{
  
  NumericMatrix accepted_samples(n, p);
  RcppParallel::RMatrix<double> w(accepted_samples);
  int count = 0;
  double dist;
  NumericVector theta_sim(p);
  
  #pragma omp parallel num_threads(ncores)
  {
    #pragma omp for
    for (int i=0; i<n; i++) {
      theta_sim[0] = unif_sitmo(i);
      theta_sim[1] = unif_sitmo(i+n);
      
      NumericVector I = SIR(762.0/763.0, 1.0/763.0, 0.0, theta_sim[0], theta_sim[1]);
      dist = calc_dist(I, x);
      if (dist < eps) {
        w(i, 0) = theta_sim[0];
        w(i, 1) = theta_sim[1];
      }
    }
  }
  
  for (int i=0; i<n; i++) {
    if (w(i, 0) != 0) {
      count += 1;
    }
  }
  std::cout << "Acceptance rate: " << (double)count / n << std::endl;
  return accepted_samples;
  
}

