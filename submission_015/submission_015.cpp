#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <sys/time.h>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <set>
#include <iomanip>

#define FOR(i,a,b)  for(__typeof(b) i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define FOREACH(x,c)   for(__typeof(c.begin()) x=c.begin();x != c.end(); x++)
#define ALL(c)      c.begin(),c.end()
#define CLEAR(c)    memset(c,0,sizeof(c))
#define SIZE(c) (int) ((c).size())

#define PB          push_back
#define MP          make_pair
#define X           first
#define Y           second

#define ULL         unsigned long long
#define LL          long long
#define LD          long double
#define II         pair<int, int>
#define DD         pair<double, double>

#define VC	    vector
#define VI          VC<int>
#define VVI         VC<VI>
#define VD          VC<double>
#define VS          VC<string>
#define VII         VC<II>
#define VDD         VC<DD>

#define DB(a)       cerr << #a << ": " << a << endl;

using namespace std;

template<class T> void print(VC < T > v) {cerr << "[";if (SIZE(v) != 0) cerr << v[0]; FOR(i, 1, SIZE(v)) cerr << "," << v[i]; cerr << "]\n";}
template<class T> string i2s(T &x) { ostringstream o; o << x; return o.str(); }
VS split(string &s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < SIZE(s)) all.PB(s.substr(p)); return all;}

//#define LOCAL

class NationSave{
    public:
        double beta, eta, alpha, delta, rho, sigma;
        int N, T;

        VD C;
        VD C0;

        VD ConsumptionDecisionRule(VD &Kt, VD &Zt){
            C0 = C;
            REP(i,N){
                //cerr << i << endl;
                double K = Kt[i];
                double Z = Zt[i];
                double W = (1-delta)*K + Z*pow(K,alpha);
                double mul1 = beta*(1-delta+Z*alpha*pow(K,alpha-1));
                double mul2 = pow(mul1,1/eta);
                C[i] = C0[i]*mul2;
                if (C[i] == 0.0 || C[i] >= 0.5*W) C[i] = W/1000;
            }
            //REP(i,N) cerr << C[i] << " "; 
            //cerr << endl;
            return C;
        }
        int SetEconomyParameters(double _beta, double _eta, double _alpha, double _delta, double _rho, double _sigma, int _N, int _T){
            beta = _beta;
            eta = _eta;
            alpha = _alpha;
            delta = _delta;
            rho = _rho;
            sigma = _sigma;
            N = _N;
            T = _T;
            C = VD(N,0.001);
            return 0;
        }
};

#ifdef LOCAL

int main(int argc, char *argv[]){
    int C;
    cin >> C;
    REP(c,C){
        double beta, eta, alpha, delta, rho, sigma;
        int N,T;
        cin >> beta >> eta >> alpha >> delta >> rho >> sigma;
        cin >> N >> T;
        cerr << beta << " " << eta << " " << alpha << " " << delta << " " << rho << " " << sigma << endl;
        cerr << N << " " << T << endl;

        NationSave ns;

        ns.SetEconomyParameters( beta, eta, alpha, delta, rho, sigma, N, T);

        VD Z(N), K(N);
        VD C;

        REP(t,T){
            REP(i,N) cin >> K[i];
            REP(i,N) cin >> Z[i];
            C = ns.ConsumptionDecisionRule(K,Z);
            REP(i,N) 
                cout << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << C[i] << endl;
            cout.flush();
        }
    }
    return 0;
}

#endif