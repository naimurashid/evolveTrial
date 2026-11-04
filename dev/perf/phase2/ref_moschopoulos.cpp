// Independent reference for gate G1 (step P2-1): the Moschopoulos (1985) series for the distribution of a sum of
// independent gammas, in long double. Not part of the package; dev/perf/phase2/g1_check.R compiles it with sourceCpp.
//
// S = sum_k theta_k G_k, G_k ~ Gamma(a_k, 1). With theta_1 = min theta, rho = sum a_k,
//   C = prod_k (theta_1 / theta_k)^a_k,  g_i = sum_k a_k (1 - theta_1 / theta_k)^i / i,
//   d_0 = 1,  d_{j+1} = 1 / (j + 1) sum_{i=1}^{j+1} i g_i d_{j+1-i},
//   Pr(S <= c) = sum_j p_j P(rho + j, c / theta_1),  p_j = C d_j  (p is a probability mass function).
// Truncating after N terms leaves sum_{j>N} p_j P(rho + j, x) <= (1 - sum_{j<=N} p_j) P(rho + N + 1, x), since the
// regularized lower gamma P(rho + j, x) decreases in j. That bound is the stopping rule.
#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <cmath>

// [[Rcpp::export]]
Rcpp::List moschopoulos_tail(Rcpp::NumericVector a, Rcpp::NumericVector theta, double c, int nmax, double tol) {
  const int K = a.size();
  int k1 = 0;
  for (int k = 1; k < K; k++) if (theta[k] < theta[k1]) k1 = k;
  const long double t1 = theta[k1];
  long double logC = 0.0L, rho = 0.0L;
  std::vector<long double> r(K);
  for (int k = 0; k < K; k++) {
    rho += a[k];
    logC += (long double) a[k] * logl(t1 / (long double) theta[k]);
    r[k] = 1.0L - t1 / (long double) theta[k];
  }
  const long double C = expl(logC);
  const double x = c / (double) t1;
  std::vector<long double> g(nmax + 2, 0.0L), d(nmax + 2, 0.0L), rp(K, 1.0L);
  d[0] = 1.0L;
  long double F = 0.0L, mass = 0.0L, bound = 1.0L;
  int j = 0;
  bool conv = false;
  for (j = 0; j <= nmax; j++) {
    if (j > 0) {
      long double gi = 0.0L;
      for (int k = 0; k < K; k++) {
        rp[k] *= r[k];
        gi += (long double) a[k] * rp[k];
      }
      g[j] = gi / j;
      long double s = 0.0L;
      for (int i = 1; i <= j; i++) s += i * g[i] * d[j - i];
      d[j] = s / j;
    }
    const long double pj = C * d[j];
    mass += pj;
    F += pj * (long double) R::pgamma(x, (double) (rho + j), 1.0, 1, 0);
    const long double rest = 1.0L - mass;
    bound = (rest > 0 ? rest : 0.0L) * (long double) R::pgamma(x, (double) (rho + j + 1), 1.0, 1, 0);
    if (bound < tol) {
      conv = true;
      break;
    }
  }
  return Rcpp::List::create(Rcpp::Named("tail") = (double) (1.0L - F), Rcpp::Named("bound") = (double) bound,
                            Rcpp::Named("terms") = j + 1, Rcpp::Named("converged") = conv);
}
