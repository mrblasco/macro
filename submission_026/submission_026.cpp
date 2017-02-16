#include <bits/stdc++.h>
using namespace std;

class NationSave
{
public:
    int n,t;
    double b,e,a,d,r,s;
    vector <double> w,cr;


    vector <double> ConsumptionDecisionRule(vector <double> Kt, vector <double> Zt)
    {
        for (int i = 0; i < n; i++) {
            w[i] = cwage(Kt[i], Zt[i]);
        }

        for (int i = 0; i < n; i++) {
            cr[i] =w[i]*d;
        }

        return cr;
    }

    int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T)
    {
        b = beta; e = eta; a = alpha;d = delta;r = rho;s = sigma;n = N;t = T;
        w.assign(N, 0.0);
        cr.assign(N, 0.0);
        return 0;
    }

    double cwage(double kt, double zt)
    {
        double ret = (1 - d) * kt + zt * pow(kt, a);
        return ret;
    }

    NationSave()
    {
    }
};


