#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <set>
#include <map>
using namespace std;


#define FOR(i, a, b)	for(int i = int(a); i < int(b); ++i)
#define FOR_(i, a, b)	for(int i = int(a); i > int(b); --i)
#define REP(i,n)	FOR(i, 0, n)
#define VC		vector
#define SZ		size()
#define DIM		resize
#define PB		push_back
#define VI		VC <int>
#define VVI		VC < VI >
#define Vc		VC <char>
#define VVc		VC < Vc >
#define VD		VC < double >
#define VVD		VC < VD >
#define PII		pair <int, int >
#define MII		map <int, int>
#define VP		VC < PII >
#define VS		VC <string>
#define VVS		VC < VS >
#define VF		VC <float>
#define VVF		VC < VF >
#define VVVF	VC < VVF >

////////////  utilities
template <class T> void _msg(string s, T x){
	cerr<<(s + " = ") << x << "\n";
}

///////////////// utilities end

class NationSave {
public:
	NationSave() {test = 0;}
	int SetEconomyParameters(
		double _beta,
		double _eta,
		double _alpha,
		double _delta,
		double _rho,
		double _sigma,
		int _N,
		int _T) {
			beta = _beta;
			eta = _eta;
			alpha = _alpha;
			delta = _delta;
			rho = _rho;
			sigma = _sigma;
			N = _N;
			T = _T;
			t = 0;
			K.clear();
			Z.clear();
			C.clear();
			W.clear();
			Coef.clear();
			Eps.clear();
			K.DIM(T + 1, VD(N));
			Z.DIM(T + 1, VD(N));
			C.DIM(T + 1, VD(N));
			W.DIM(T + 1, VD(N));
			Coef.DIM(T + 1, VD(N));
			Eps.DIM(T + 1, VD(N));
			eps_t.clear();
			eps_t.DIM(T + 1);
			eps_tot = 0;
			test++;
			modif = 0;
			ups = 0;
			return 0;
		}
	VD ConsumptionDecisionRule(VD Kt, VD Zt) {
		REP(n, N) {
			K[t][n] = Kt[n];
			Z[t][n] = Zt[n];
			W[t][n] = (1 - delta) * K[t][n] + Z[t][n] * pow(K[t][n], alpha);
			if(t == 0) {
				C[t][n] = 0.001 * W[t][n];
			} else if(t < T - 1) {
			    Coef[t][n] = coef(n);
				C[t][n] = Coef[t][n] * C[t - 1][n];
				if(C[t][n] >= W[t][n]) {
					++modif;	//
					if(delta < 0.5)
						C[t][n] = W[t][n] * (delta);
					else
						C[t][n] = W[t][n] * (1 - delta);
				}
				if(C[t][n] < 1e-300) {
					C[t][n] = 0.001 * W[t][n];
					ups ++;	//
				}
				Eps[t][n] = scor(n);
				eps_t[t] += Eps[t][n];
			} else {	// t = T - 1
				Coef[t][n] = coef(n);
				if(n < 0.6 * N) {
					C[t][n] = 0.7 * W[t][n];
				} else {
					C[t][n] = Coef[t][n] * C[t - 1][n];
					if(C[t][n] >= W[t][n]) {
						++modif;	//
						C[t][n] = 0.7 * W[t][n];
					}
					if(C[t][n] < 1e-300) {
						C[t][n] = 0.7 * W[t][n];
						ups ++;	//
					}
					
				}
				Eps[t][n] = scor(n);
				eps_t[t] += Eps[t][n];
			}
		}
		if(t) {
			compensate();
			eps_t[t] = 0;
			REP(n, N) {
				Eps[t][n] = scor(n);
				eps_t[t] += Eps[t][n];
			}
			eps_tot += fabs(eps_t[t] / N);
		}
		res.clear();
		REP(n, N)
			res.PB(C[t][n]);
		++t;
		/*
		if(t == T) {
			_msg("modif", modif);
			_msg("ups", ups);
			_msg("eps_tot", eps_tot);
			double eps_max = 0;
			REP(ii, T) {
                if(fabs(eps_t[ii]) > eps_max)
                    eps_max = fabs(eps_t[ii]);
			}
			_msg("eps_max", eps_max);
		}
		*/
		return res;
	}
private:
	double scor(int i) {
		return beta * pow(C[t - 1][i] / C[t][i], eta) * (1 - delta + Z[t][i] * alpha * pow(K[t][i], alpha - 1)) - 1;
	}
	double coef(int i) {
		return pow((beta * (1 - delta + Z[t][i] * alpha * pow(K[t][i], alpha - 1))), 1.0 / eta);
	}
	void compensate() {
	    VI mark(N);
	    VD E1 = Eps[t];
	    VD C1 = C[t];
	    int poz = 0;
	    REP(i, N)
            if(E1[i] >= 0) {
                mark[i] = 1;
                poz++;
            }
        if(poz == 0) // poz < 50
            return;
	    double eps_med = eps_t[t] / poz;	// negative
	    REP(i, N)
            if(mark[i]) {
                E1[i] = E1[i] - eps_med;
                C1[i] = C_calc(i, E1[i]);
            }
		REP(i, N)
            if(mark[i]) {
				C[t][i] = C1[i];
				Eps[t][i] = E1[i];
			}
	}
	double C_calc(int i, double e) {
		return C[t - 1][i] * pow(beta * (1 - delta + Z[t][i] * alpha * pow(K[t][i], alpha - 1)) / (1 + e), 1.0 / eta);
	}
	
	double beta;
	double eta;
	double alpha;
	double delta;
	double rho;
	double sigma;
	int N;
	int T;
	VD res;
	int t;
	int test;
	VVD K, Z, C, W;
	int modif;
	int ups;
	VVD Coef, Eps;
	VD eps_t;
	double eps_tot;
};
// moshu