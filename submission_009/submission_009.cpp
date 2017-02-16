#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>


static int seed = -1;

#define MAX_TIME 120

#ifdef DEBUG
#define RATIO 0.4
#else
#define RATIO 1.0
#endif


#define M 10000

#define SIGMA 0.01

double AMP = 0.01;


#define ALL(x) std::begin(x), std::end(x)
#define X(x) ((x) / N_)
#define Y(x) ((x) % N_)
#define XY(x, y) (N_ * (x) + (y))


class xor128_t {
public:
  xor128_t() : s0_(1), s1_(2) {};

public:
  inline uint64_t operator()() {
    uint64_t s1 = s0_;
    uint64_t s0 = s0_ = s1_;

    s1 ^= s1 << 23;
    s1 ^= s1 >> 17;
    s1 ^= s0;
    s1 ^= s0 >> 26;

    return s0_ + (s1_ = s1);
  };

private:
  uint64_t s0_;
  uint64_t s1_;
};


xor128_t rng;


class rnorm_t {
public:
  rnorm_t() : c_(0) {};
  
public:
  template<class T>
  inline double operator()(T& rng) {
    if ((c_ = 1 - c_))
      return X_;

    double U1, U2, W, x;

    do {
      U1 = 2. * rng() / UINT64_MAX - 1.0;
      U2 = 2. * rng() / UINT64_MAX - 1.0;

      W = U1 * U1 + U2 * U2;
    } while (! (0.0 < W && W < 1.0));

    X_ = U2 * (x = sqrt((-2.0 * log(W)) / W));

    return U1 * x;
  };

private:
  int c_;

  double X_;
};


rnorm_t rnorm;


inline unsigned long long rdtsc()
{
#ifdef __amd64
  unsigned long long a, d;

  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));

  return (d << 32) | a;
#else
  unsigned long long x;

  __asm__ volatile ("rdtsc" : "=A" (x));

  return x;
#endif
}

inline double gettime()
{
  return rdtsc() / 2500193521.0;
}


class z_sampler_t {
public:
  z_sampler_t(double rho = 1.0, double sigma = 1.0) : rho_(rho), sigma_(sigma) {};

public:
  double operator()(const double z, int N = 32) const {
    double s = 0.0, rhologz(rho_ * log(z));

    for (int i = 0; i < N; i ++)
      s += exp(rhologz + rnorm(rng) * sigma_);

    return s / N;
  };

private:
  double rho_;
  double sigma_;
};


class NationSave {
public:
  NationSave() : st_(gettime()), t_(0), quiet_(false) {
    std::cerr << std::fixed << std::setprecision(2);
  };

  ~NationSave() {
    delete[] K_;
    delete[] Z_;
    delete[] W_;
    delete[] C_;
  };

public:
  int SetEconomyParameters(double beta,
                           double eta,
                           double alpha,
                           double delta,
                           double rho,
                           double sigma,
                           int N,
                           int T) {

    beta_  = beta;
    eta_   = eta;
    alpha_ = alpha;
    delta_ = delta;
    rho_   = rho;
    sigma_ = sigma;
    
    N_ = N;
    T_ = T;

    std::cerr << "solver: seed=" << seed << " " <<
      "beta="  << beta_  << ' ' <<
      "eta="   << eta_   << ' ' <<
      "alpha=" << alpha_ << ' ' <<
      "delta=" << delta_ << ' ' <<
      "rho="   << rho_   << ' ' <<
      "sigma=" << sigma_ << ' ' <<
      "N="     << N_     << ' ' <<
      "T="     << T_     << std::endl;


    int size = (T_ + 1) * N_;
    
    K_ = new double[size];
    Z_ = new double[size];
    W_ = new double[size];
    C_ = new double[size];

    std::fill(K_, K_ + size, 1.0);
    std::fill(Z_, Z_ + size, 1.0);

#ifdef DEBUG
    {
      double *p0 = K_, *p1 = Z_, *p2 = W_, *p3 = C_;
    
      for (int t = 0; t <= T_; t ++)
        for (int i = 0; i < N_; i ++) {
          assert(*p0 ++ == 1.0);
          assert(*p1 ++ == 1.0);
          assert(*p2 ++ == 0.0);
          assert(*p3 ++ == 0.0);
        }
    }
#endif

    z_sample_ = z_sampler_t(rho_, sigma_);

    T_ = std::max(T_ / 10, 200);

    quiet_ = true;

    double sp = 0.0, ampp = -1.0, s;

    for (int x = 0; AMP = x / 100.0 + 0.01, x < 60; x += 2) {
      std::vector<double> Kt(N_, 1.0), Zt(N_, 1.0), Wt(N_);

      for (int t = 0; t < T_; t ++) {
        if (t > 0)
          for (int i = 0; i < N_; i ++)
            Zt[i] = exp(rho_ * log(Zt[i]) + rnorm(rng) * sigma_);

        for (int i = 0; i < N_; i ++)
          Wt[i] = (1.0 - delta_) * Kt[i] + Zt[i] * pow(Kt[i], alpha_);
        
        it_ = gettime();
        
        std::vector<double> Ct(ConsumptionDecisionRule(Kt, Zt));

        for (int i = 0; i < N_; i ++) {
          assert(Ct[i] < Wt[i]);
          
          Kt[i] = Wt[i] - Ct[i];
        }
      }

      if ((s = done()) > sp) {
        std::cerr << "pretest: seed=" << seed << ' ' <<
          "x="   << x   << ' ' <<
          "s="   << s   << ' ' <<
          "amp=" << AMP << '*' << std::endl;
        
        sp   = s;
        ampp = AMP;
      }

      t_ = 0;

      if (s < 970000.0)
        break;
    }

    AMP = ampp;

    T_ = T;

    quiet_ = false;

    it_ = gettime();

    return 0;
  };
 
  std::vector<double> ConsumptionDecisionRule(std::vector<double> Kt,
                                              std::vector<double> Zt) {

    st_ += gettime() - it_;

    if (! quiet_)
      if (t_ % 100 == 0)
        std::cerr << "solver: seed=" << seed << ' ' <<
          "t="       << t_                        << ' ' <<
          "elapsed=" << (gettime() - st_) / RATIO << std::endl;


    assert(Kt.size() == N_);
    assert(Zt.size() == N_);

    for (int i = 0, p = XY(t_, 0); i < N_; i ++, p ++) {
      assert(fabs(Kt[i] - K_[p]) < 1.0e-6);
      
      K_[p] = Kt[i];
      Z_[p] = Zt[i];
    }
    

    double *Cp = new double[N_], *zpp = new double[N_];
    
    
    for (int i = 0, p = XY(t_, 0); i < N_; i ++, p ++)
      W_[p] = (1.0 - delta_) * K_[p] + Z_[p] * pow(K_[p], alpha_);

    for (int i = 0, p = XY(t_, 0); i < N_; i ++, p ++)
      Cp[i] = C_[p] = W_[p] * AMP;

    for (int i = 0, p = XY(t_, 0), pp = p + N_; i < N_; i ++, p ++, pp ++) {
      K_[pp] = W_[p] - C_[p];

      zpp[i] = Z_[pp] = z_sample_(Z_[p]);
        
      W_[pp] = (1.0 - delta_) * K_[pp] + Z_[pp] * pow(K_[pp], alpha_);

      if (t_ < T_ - 1) {
        C_[pp] = W_[pp] * AMP;
      }
      else {
        C_[pp] = W_[pp];
      }
    }
      
    double sp = 1.0e+30, s0 = 0.0, s1 = 0.0, e0[222] = {0}, e1[222];

    if (t_ > 0)
      for (int i = 0; i < N_; i ++)
        s0 += e0[i] = eer_t_i(t_ - 1, i);
      
    for (int i = 0; i < N_; i ++)
      s1 += e1[i] = eer_t_i(t_, i);

    for (int m = 0, MM = t_ < T_ - 1 ? M : M * 10; m < MM; m ++) {
      int I = rng() % N_, p = XY(t_, I), pp = p + N_;
        
      double a[4];

      a[0] = C_[p ];
      a[1] = K_[pp];
      a[2] = W_[pp];
      a[3] = C_[pp];

      double c;

      do {
        c = C_[p] + W_[p] * rnorm(rng) * SIGMA;
      } while (! (0.0 < c && c < W_[p]));

      K_[pp] = W_[p] - (C_[p] = c);

      W_[pp] = (1.0 - delta_) * K_[pp] + Z_[pp] * pow(K_[pp], alpha_);

      if (t_ < T_ - 1) {
        C_[pp] = W_[pp] * AMP;
      }
      else {
        C_[pp] = W_[pp];
      }

      double ss0 = 0.0, ss1 = 0.0, ee0 = 0.0, ee1 = 0.0;

      if (t_ > 0)
        ss0 = s0 - e0[I] + (ee0 = eer_t_i(t_ - 1, I));

      ss1 = s1 - e1[I] + (ee1 = eer_t_i(t_, I));
                            
      double s = fabs(ss0) + fabs(ss1);

      if (s < sp) {
#ifdef DEBUG
#if 0
        std::cerr << "solver: seed=" << seed << ' ' <<
          "t="   << t_                                                  << ' ' <<
          "m="   << m                                                   << ' ' <<
          "I="   << I                                                   << ' ' <<
          "eer=" << s / N_ << '(' << ss0 / N_ << '/' << ss1 / N_ << ')' << '*' << std::endl;
#endif
#endif
            
        sp = s;

        s0 = ss0;
        s1 = ss1;

        e0[I] = ee0;
        e1[I] = ee1;

        memcpy(Cp, C_ + XY(t_, 0), sizeof(double) * N_);
      }
      else {
        C_[p ] = a[0];
        K_[pp] = a[1];
        W_[pp] = a[2];
        C_[pp] = a[3];
      }
    }

    for (int i = 0, p = XY(t_, 0), pp = p + N_; i < N_; i ++, p ++, pp ++)
      K_[pp] = W_[p] - (C_[p] = Cp[i]);

    
    std::vector<double> a(N_);

    for (int i = 0, p = XY(t_, 0); i < N_; i ++, p ++)
      a[i] = C_[p];


    delete[] Cp;
    delete[] zpp;

    t_ ++;

    it_ = gettime();

    return std::move(a);
  };

public:
  double done() {
    assert(t_ == T_);

    for (int i = 0, pp = XY(T_, 0), p = pp - N_; i < N_; i ++, p ++, pp ++) {
      Z_[pp] = z_sample_(Z_[p]);

      C_[pp] = W_[pp] = (1.0 - delta_) * K_[pp] + Z_[pp] * pow(K_[pp], alpha_);
    }

#ifdef DEBUG
#if 0
    std::cerr << std::setprecision(6);
    
    for (int t = 0; t < T_; t ++)
      std::cerr << "solver: seed=" << seed << ' ' <<
        "t="     << t                                 << ' ' <<
        "eer_t=" << eer_t(t) << std::endl;

    std::cerr << std::setprecision(2);
#endif
#endif

    double s = eer(), score = 1.0 / (1.0 + s), score2 = 1.0 - sqrt(1.0 - score);

    score  *= 1.0e+6;
    score2 *= 1.0e+6;

    if (! quiet_)
      std::cerr << "solver: seed=" << seed << ' ' <<
        "eer="   << std::setprecision(6) << s     << ' ' <<
        "amp="   << std::setprecision(2) << AMP   << ' ' <<
        "score=" << std::setprecision(2) << score << '(' << score2 << ')' << std::endl;

#if 0
    for (int t = 0; t < T_; t ++) {
      std::cerr << "t=" << t << ": ";
      for (int i = 0, p = XY(t, i), pp = p + N_; i < 8; i ++, p ++, pp ++)
        std::cerr <<
          C_[p] << ' ' << C_[pp] << ' ' << Z_[pp] << ' ' << K_[pp] << "    ";
      std::cerr << eer_t(t) << std::endl;
    }
#endif

    return score;
  };

private:
  double eer_t_i(int t, int i) const {
    const int p = XY(t, i), pp = p + N_;
  
    double U = pow(C_[p] / C_[pp], eta_);
      
    double V = 1.0 - delta_;

    double W = Z_[pp] * alpha_ * pow(K_[pp], alpha_ - 1.0);

    assert(U >= 0.0);
    assert(V >= 0.0);
    assert(W >= 0.0);

    assert(beta_ * U * (V + W) >= 0.0);

    return beta_ * U * (V + W) - 1.0;
  };

  double eer_t(int t) const {
    double ss = 0.0;

    for (int i = 0; i < N_; i ++)
      ss += eer_t_i(t, i);

    return ss / N_;
  };

  double eer() {
    const int SAMPLES = 2000;
  
    double* stat = new double[(T_ + 1) * N_];

    double* table = new double[T_];

    double s = 0.0;

    for (int t = 0; t < T_; t ++) {
      double ss = 0.0;

      for (int i = 0; i < N_; i ++)
        ss += stat[i * T_ + t] = eer_t_i(t, i);

      s += table[t] = fabs(ss / N_);
    }

    if (T_ < SAMPLES) {
      s = table[T_ - 1];
    
      for (int i = 0; i < SAMPLES - 1; i ++)
        s += table[rng() % (T_ - 1)];
    }

#if 0
    s = 0.0;
    
    for (int i = 0; i < T_ - 1; i ++)
      s += table[i];

    s /= (T_ - 1);
#endif

    for (int i = 0; i < N_; i ++) {
      double *p = stat + i * T_,  mean = std::accumulate(p, p + T_, 0.0) / T_, var = 0.0, d;

      for (int t = 0; t < T_; t ++, p ++) {
        d = *p - mean;
      
        var += d * d;
      }

      var /= T_;

#if 0
      std::cerr << "score: seed=" << seed << ' ' <<
        "i="      << i         << ' ' <<
        "mean*="  << mean      << ' ' <<
        "sigma*=" << sqrt(var) << std::endl;
#endif
    }

    delete[] stat;
    delete[] table;

    return s / SAMPLES;
  };

private:
  double beta_;
  double eta_;
  double alpha_;
  double delta_;
  double rho_;
  double sigma_;
  
  int N_;
  int T_;

private:
  int t_;

private:
  double* K_;
  double* Z_;
  double* W_;
  double* C_;

private:
  z_sampler_t z_sample_;

private:
  double st_;
  double it_;

private:
  bool quiet_;
};
