#include <omp.h>
#include <sitmo.h>
#include <cmath>
using namespace Rcpp;


// Function to simulate data from SIR model

double unif_sitmo(int seed) {
  uint32_t coreseed = static_cast<uint32_t>(seed);
  sitmo::prng eng(coreseed);
  double mx = sitmo::prng::max();
  double x = eng() / mx;
  return x;
}

double calc_dist(double* x_sim, NumericVector x) {
  double total;
  #pragma omp parallel reduction( + : total)
  {
    double inner_sum;
    #pragma omp for
    for (int i=0; i<20; i++) {
      inner_sum += pow(x[i]-x_sim[i], 2);
    }
    total += inner_sum;
  }
  return total;
}

double calc_dist_serial(double* x_sim, NumericVector x) {
  double total;
  for (int i=0; i<20; i++) {
    total += pow(x[i]-x_sim[i], 2);
  }
  return total;
}

double ABC(int n, double eps, int p)
{
  int count = 0;
  double dist;

  #pragma omp parallel num_threads(8)
  {
    double theta_sim[2];
    #pragma omp for
    for (int i=0; i<n; i++) {

      #pragma omp critical
      {
      theta_sim[0] = unif_sitmo(i);
      theta_sim[1] = unif_sitmo(i+n);
      }
      
      // NumericVector I = SIR(762.0/763.0, 1.0/763.0, 0.0, theta_sim[1], theta_sim[2]);

      double S[20];
      double I[20];
      double R[20];

      S[0] = 762.0/763.0;
      I[0] = 1.0/763.0;
      R[0] = 0.0;

      for (int t = 1; t < 20; t++) {
        S[t] = S[t-1] - theta_sim[0] * S[t-1] * I[t-1];
        I[t] = I[t-1] + theta_sim[0] * S[t-1] * I[t-1] - theta_sim[1] * I[t-1];
        R[t] = R[t-1] + theta_sim[1] * I[t-1];
      }

      #pragma omp critical
      {
      dist = calc_dist_serial(I, x);
      if (dist < eps) {
        count++;
      }
      }
    }
  }
  
  double acceptance_rate = (double) count / n;
  std::cout << eps << " " << acceptance_rate << std::endl;
  return acceptance_rate;
  
}

int main() {
  double eps_min = 0;
  double eps_max = 5;
  double eps_step = 0.05;

  double acceptance_rates[(eps_max-eps_min)/eps_step];
  double ar;

  int n = 100;
  for (double eps = eps_min; eps < eps_max; eps += eps_step) {
    ar = ABC(n, eps, 2);
    acceptance_rates[(eps-eps_min)/eps_step] = ar;
  }

  return 0;
}