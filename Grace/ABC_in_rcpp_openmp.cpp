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
NumericVector SIR(NumericVector s, NumericVector i, NumericVector r, double s0, double i0, double r0, double beta, double gamma) {
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

// Function to calculate distance
// [[Rcpp::export]]
double calc_dist(NumericVector x_sim, NumericVector x) {
  if (x_sim.size() != x.size()) {
    stop("Vectors must be of the same length");
  }
  
  double sum = 0.0;
  int n = x.size();

  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; ++i) {
    sum += std::pow(x_sim[i] - x[i], 2);
  }
  
  return std::sqrt(sum);
}

// [[Rcpp::export]]
NumericMatrix ABC(int n, double eps, int p, NumericVector x, int ncores)
{
  
  NumericMatrix accepted_samples(n, p);
  int count = 0;
  NumericVector zeros(x.size(), 0);

  #pragma omp parallel num_threads(ncores)
  {
    NumericVector theta_sim(p);
    NumericVector x_sim(x.size());
    int local_count = 0;
    
    #pragma omp for
    for (int i = 0; i < n; i++) {
      theta_sim[0] = unif_sitmo(i); // sample beta from prior
      theta_sim[1] = unif_sitmo(i + n); // sample gamma from prior
      x_sim = SIR(zeros, zeros, zeros, 762.0/763.0, 1.0/763.0, 0, theta_sim[0], theta_sim[1]); // simulate data from SIR model
      #pragma omp critical
      {
      std::cout << "zeros: " << zeros << std::endl;
      std::cout << "theta_sim[0]: " << theta_sim[0] << std::endl;
      std::cout << "theta_sum[1]: " << theta_sim[1] << std::endl;
      std::cout << "x_sim: " << x_sim << std::endl;
      }
      double dist = calc_dist(x_sim, x); // calculate distance
      //std::cout << "dist: " << dist << std::endl;
      // Accept-reject step
      if (dist <= eps) {
        #pragma omp critical
        {
          for (int j = 0; j < p; j++) {
            accepted_samples(count, j) = theta_sim[j];
          }
          count++;
        }
      }
    }
  }
  
  Rprintf("Acceptance rate: %f\n", (double)count / n);
  return accepted_samples;
  
}
