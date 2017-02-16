#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <tr1/random>
#include <sys/time.h>
#include <float.h>
#include <iomanip>
#define SETCOL  setw(20) << left
using namespace std;



class NationSave {
	private:
	double beta, eta, alpha, delta, rho, sigma;
	int N , T, iteration;
	int start;
	vector< vector<double> > Kti;
	vector< vector<double> > Zti;
	vector< vector<double> > Cti;
	vector< vector<double> > eerti;
	vector<double> Cprev;

	int getTime() {
		timeval tv;
		gettimeofday(&tv, 0);
		return 1000000*tv.tv_sec + tv.tv_usec;
	}

	double calcWage(double K, double Z) {
		return (1 - delta)*K + Z*pow(K, alpha);
	}

	bool simulate(int simulation, double f, double Wage) {
		Cti[iteration][simulation] = f*Wage;
		for(int i=iteration; i<=min(T, iteration+T/10); i++) {

			if(i > iteration) Zti[i][simulation] = exp(rho * log(Zti[i - 1][simulation]));

			Wage = calcWage(Kti[i][simulation], Zti[i][simulation]);

			if(i == iteration) Cti[i][simulation] = f*Wage;

			if(i > iteration) Cti[i][simulation] = Cti[i-1][simulation] / pow( 1/(beta * ( 1 - delta + Zti[i][simulation]*alpha*pow(Kti[i][simulation], alpha - 1) )), 1/eta );
			//cerr << simulation << endl;
			if(Cti[i][simulation] > Wage) return false;

			if(i+1 <= T) Kti[i+1][simulation] = Wage - Cti[i][simulation];
			//if(i == iteration) cerr << "f: " << SETCOL << f << "sim: " << SETCOL << simulation << "K: " << SETCOL << Kti[i][simulation] << "Z: " << SETCOL << Zti[i][simulation] << "W: " << SETCOL << Wage << "C: " << SETCOL << Cti[i][simulation] << endl;

		}
		return true;
	}
	double calcscoreoveronesimulation(int simulation) {

		for(int t=0; t<T; t++) {
			eerti[t][simulation] = beta*( pow(Cti[t][simulation]/Cti[t+1][simulation], eta) * ( 1 - delta + Zti[t+1][simulation]*alpha*pow(Kti[t+1][simulation], alpha - 1) ) ) - 1;
		}

		double eer = 0;
		for(int i=0; i<iteration-1; i++) {
			eer += eerti[i][simulation];
		}
		//cerr << eer << endl;

		return abs(eer);
	}


	public:

	int SetEconomyParameters(double Beta, double Eta, double Alpha, double Delta, double Rho, double Sigma, int NN, int TT) {
		start = getTime();

		beta = Beta;
		eta = Eta;
		alpha = Alpha;
		delta = Delta;
		rho = Rho;
		sigma = Sigma;

		cerr << "beta: " << beta << endl;
		cerr << "eta: " << eta << endl;
		cerr << "alpha: " << alpha << endl;
		cerr << "delta: " << delta << endl;
		cerr << "rho: " << rho << endl;
		cerr << "sigma: " << sigma << endl;

		N = NN;
		T = TT;
		iteration = 0;

		Kti.resize(T+1);
		Zti.resize(T+1);
		Cti.resize(T+1);
		eerti.resize(T+1);
		for(int i=0; i<T+1; i++) {
			Kti[i].resize(N);
			Zti[i].resize(N);
			Cti[i].resize(N);
			eerti[i].resize(N);
		}
		cerr.precision(10);

		return 0;
	}

	vector<double> ConsumptionDecisionRule(vector<double> Kt, vector<double> Zt) {
		vector<double> result(N);
		vector<double> Wages(N);
		Kti[iteration] = Kt; 
		Zti[iteration] = Zt; 

		if(iteration == 0) {
			double Wage = calcWage(Kt[0], Zt[0]);
			double bestf = 0.5;
/*			double f = 0.0000000001;
			double fff = 0.9999999999;
			double ff = (f + fff)/2;
			while(abs(fff - f) >= 0.0000000000001) {
				bool sim1 = simulate(0, f, Wage);
				bool sim2 = simulate(0, ff, Wage);
				bool sim3 = simulate(0, fff, Wage);
				if(sim1 == true && sim2 == true && sim3 == false) {
					bestf = ff;
					f = ff;
					ff = (f + fff)/2;
				}
				if(sim1 == true && sim2 == false && sim3 == false) {
					bestf = f;
					fff = ff;
					ff = (f + fff)/2;
				}
				if(sim1 == true && sim2 == true && sim3 == true) {
					bestf = f;
					break;
				}
			}*/
			bestf = 0.1;

			//cerr << SETCOL << bestf << endl;

			for(int i=0; i<N; i++) {
				double C = bestf * Wage;
				result[i] = C;
			}
		} else {

			for(int i=0; i<N; i++) {
				double Wage = calcWage(Kt[i], Zt[i]);
				Wages[i] = Wage;

				double C = Cprev[i] / pow( 1/(beta * ( 1 - delta + Zt[i]*alpha*pow(Kt[i], alpha - 1) )), 1/eta );
				//if(C == 0) C = DBL_MIN;

				if(C >= Wage || C < 2*DBL_MIN) {
					double bestf = 0.5;
					double f = 0.0000000001;
					double fff = 0.9999999999;
					double ff = (f + fff)/2;
					while(abs(fff - f) >= 0.0000000001) {
						bool sim1 = simulate(i, f, Wage);
						bool sim2 = simulate(i, ff, Wage);
						bool sim3 = simulate(i, fff, Wage);
						if(sim1 == true && sim2 == true && sim3 == false) {
							bestf = ff;
							f = ff;
							ff = (f + fff)/2;
						}
						if(sim1 == true && sim2 == false && sim3 == false) {
							bestf = f;
							fff = ff;
							ff = (f + fff)/2;
						}
						if(sim1 == true && sim2 == true && sim3 == true) {
							bestf = f;
							break;
						}
					}
					//cerr << SETCOL << bestf << SETCOL << iteration << i << endl;
					C = bestf*Wage;
				}
				result[i] = C;
			}
		}

		Cti[iteration] = result;
		Cprev = result;

		iteration++;
		if(iteration > 1999) cerr << (getTime() - start)/1000000.0 << endl;
		return result;
	}

};

