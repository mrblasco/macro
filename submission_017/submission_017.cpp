#include<cstdio>
#include<iostream>
#include<vector>
#include<cmath>

using namespace std;



double calculateErr(double beta, double eta, double delta, double alpha, double Kt, double Zt, double wage, double Ct, double Ct1){
    double err;
    err = beta * (pow(Ct/Ct1, eta) * (1.0 - delta + Zt * alpha * pow(Kt, alpha - 1.0))) -1;
    if(err < 0){
        err = -err;
    }
    return err;
}

double calC(double beta, double eta, double delta, double alpha, double Kt, double Zt, double Ct, double wage, double kstar){
    double y = 1.0 / beta;
    double d = 1.0 -delta + Zt * alpha * pow(Kt, alpha-1.0);
    y = y / d;
    double res, tmp;
    tmp = pow(y, 1.0 / eta);
    if(tmp == 0){
        res = 0.37 * wage;
        res = wage - kstar;
    }
    else{
        res = Ct / tmp;
    }
    if(res > wage || res <= 0){
        res = 0.37 * wage;
    }
    return res;
}

double first(double beta, double eta, double alpha, double rho, double Kt, double Zt, double wage){
    double C;
    C = 0.22 * wage;
    return C;
}


class NationSave{
    public:
    double beta, eta, alpha, delta, rho, sigma;
    bool flag;
    int N, T;
    vector<double> current;
    public:
    vector<double> ConsumptionDecisionRule(vector<double> Kt, vector<double> Zt){
        vector<double> Ct;
        double tmp, wage, kstar;
        for(int i = 0; i < Kt.size(); i++){
            //Ct.push_back(0.8 * Kt[i]* (1-delta));
            wage = (1.0 - delta) * Kt[i] + Zt[i] * pow(Kt[i], alpha);
            kstar = pow(Zt[i] *alpha *  beta / (1.0 - beta + delta * beta), 1.0 / (1.0 - alpha));
            if(flag){
                tmp = 0.22 * wage;
            }
            else{
                tmp = 0.25 * wage;
                tmp = calC(beta, eta, delta, alpha, Kt[i], Zt[i], current[i], wage, kstar);
            }
            Ct.push_back(tmp);
            //cout<<flag<<"==="<<endl;
        }
        flag = false;
        current.clear();
        for(int i = 0; i < Kt.size(); i++){
            current.push_back(Ct[i]);
        }
        return Ct;
    }

    int SetEconomyParameters(double beta2, double eta2,double alpha2,double delta2,double rho2,double sigma2,int N2,int T2){
        beta = beta2;
        eta = eta2;
        alpha = alpha2;
        delta = delta2;
        rho = rho2;
        sigma = sigma2;
        N = N2;
        T = T2;
        flag = true;
        return 0;
    }
};