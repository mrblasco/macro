using System;
using System.Linq;
using System.Collections.Generic;

public class NationSave
{
    int n,t;
    double b,e,a,d,r,s;
    double[] w,cr;


    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
    {
        for (int i = 0; i < n; i++) {
            w[i] = cwage(Kt[i], Zt[i]);
        }

        for (int i = 0; i < n; i++) {
            cr[i] =w[i]*d;
        }

        return cr;
    }

    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T)
    {
        b = beta; e = eta; a = alpha;d = delta;r = rho;s = sigma;n = N;t = T;
        w = new double[N];
        cr = new double[N];
        return 0;
    }

    double cwage(double kt, double zt)
    {
        double ret = (1 - d) * kt + zt * Math.Pow(kt, a);
        return ret;
    }

    public NationSave()
    {
    }
}


