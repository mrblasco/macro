#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class NationSave
{
public:
	int SetEconomyParameters(double beta, double eta, double alfa, double delta, double rho, double sigma, int N, int T);
	vector <double> ConsumptionDecisionRule(vector<double> Kt, vector<double> Zt);

private:
	double beta, eta, alfa, delta, rho, sigma;
    double Cmin;
	int N, T;
	int num;
	vector <double> wage;
    vector <double> C;


	void calculateWage(vector <double> Kt, vector <double> Zt);

};







int NationSave::SetEconomyParameters(double beta, double eta, double alfa, double delta, double rho, double sigma, int N, int T)
{
	this->beta = beta;
	this->eta = eta;
	this->alfa = alfa;
	this->delta = delta;
	this->rho = rho;
	this->sigma = sigma;
	this->N = N;
	this->T = T;
	num = 0;
  
    Cmin=1E-300;
  	C.clear();
	return 1;
}

vector <double> NationSave::ConsumptionDecisionRule(vector<double> Kt, vector<double> Zt)
{
	if (num == 0)
	{
		for (int n = 0; n != N; n++)
			wage.push_back(0.0);
		calculateWage(Kt, Zt);
		for (int n = 0; n != N; n++)
			C.push_back(wage[n] *0.4);
        num=1;
		return C;
	}
	
	double tmp, tmpbeta, tmpinv,powtmp;
	const double inveta = 1. / eta;
	calculateWage(Kt, Zt);
	for (int n = 0; n != N; n++)
	{
      	tmp = 1. - delta + Zt[n] * alfa*pow(Kt[n], (alfa - 1.));
		tmpbeta = beta*tmp;;
		powtmp = pow(tmpbeta, inveta);
		C[n] *= powtmp;
        if(C[n]<Cmin) C[n]=Cmin;
		if (C[n] > 0.8*wage[n]) C[n] = 0.8*wage[n]; 
		if (C[n] < 0) C[n] = 0;	
      		 
	}
	return C;
	
}

void NationSave::calculateWage(vector <double> Kt, vector <double> Zt)
{

	for (int n = 0; n != N; n++)
		wage[n] = (1 - delta) * Kt[n] + Zt[n] * pow(Kt[n], alfa);
}



