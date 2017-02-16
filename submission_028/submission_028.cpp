#include <cmath>
#include <vector>
using namespace std;

// Dead simple strategy, just to see how the process of
// submitting a solution works.
class NationSave {
  double beta, eta, alpha, delta, rho, sigma;
  int N, T;

  double solve(double Kt, double Zt) {
    auto Wt = (1 - delta) * Kt + Zt * std::pow(Kt, alpha);
    return Wt / 2;
  }

public:

  int SetEconomyParameters(
      double beta, double eta, double alpha, double delta,
      double rho, double sigma, int N, int T) {
    this->beta = beta;
    this->eta = eta;
    this->alpha = alpha;
    this->delta = delta;
    this->rho = rho;
    this->sigma = sigma;
    this->N = N;
    this->T = T;
    return 0;
  }

  std::vector<double> ConsumptionDecisionRule(
      const std::vector<double>& Kt,
      const std::vector<double>& Zt) {
    std::vector<double> ret(N, 0.0);
    for (int i = 0; i < N; i++) {
      ret[i] = solve(Kt[i], Zt[i]);
    }
    return ret;
  }
};

