#include <Rcpp.h>
#include <RcppParallel.h>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]

// Function to simulate data from SIR model
// [[Rcpp::export]]
NumericVector SIR(double s0, double i0, double r0, double beta, double gamma, int T) {
  double s[T];
  double i[T];
  double r[T];
  
  s[0] = s0;
  i[0] = i0;
  r[0] = r0;

  for (int t = 1; t < T; t++) {
    s[t] = s[t-1] - beta * s[t-1] * i[t-1];
    i[t] = i[t-1] + beta * s[t-1] * i[t-1] - gamma * i[t-1];
    r[t] = r[t-1] + gamma * i[t-1];
  }
  
  NumericVector i_return;
  
  for (int t = 1; t < T; t++) {
    i_return[t] = i[t];
  }

  return i_return;
}


// Function to calculate distance
double calc_dist(NumericVector x_sim, NumericVector x) {
  double total = 0;
  for (int i=0; i<x.size(); i++) {
    total += pow(x[i]-x_sim[i], 2);
  }
  return total;
}

// [[Rcpp::export]]
NumericMatrix ABC(int n, double eps, int p, NumericVector x, int ncores) {
  
  NumericMatrix accepted_samples;
  int count = 0;

  #pragma omp parallel num_threads(ncores)
  {
  
    NumericVector theta_sim(p);
    NumericVector x_sim(p);
    double dist;

    #pragma omp for
    for(int i = 0; i < n; i++) {
      theta_sim = runif(p, 0, 1); // sample from prior
      
      x_sim = SIR(762/763, 1/763, 0, theta_sim[1], theta_sim[2], 20); // simulate data from SIR model
      dist = calc_dist(x_sim, x); // calculate distance
    
      // Accept-reject step
      if(dist <= eps){
        for(int j = 0; j < p; j++) {
          accepted_samples(i, j) = theta_sim[j];
        }
        count++;
      }
    }

  }
  
  Rprintf("Acceptance rate: %f\n", (double)count / n);
  return accepted_samples;
  
}
