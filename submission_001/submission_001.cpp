#include <bits/stdc++.h>

using namespace std;

#define PI 3.14159265
#define USE_IMPROVED_VERSION 1
#define logout cerr

class EconomyParams
{
public:
    double beta;
    double eta;
    double alpha;
    double delta;
    double rho;
    double sigma;
    int seed;
    int case_index;
    double score;
    int N;
    int T;
};

class NationSave
{
public:
    int SetEconomyParameters(double beta, double eta, double alpha,
                             double delta, double rho, double sigma,
                             int N, int T);
    vector <double> ConsumptionDecisionRule(vector <double> & Kt,
                                            vector <double> & Zt);
    double find_zero1(EconomyParams & ecoParams, double WTm1, double ZT,
                      double x0, double x1, double tol);
    double getNormalDist(const double & mu, const double & sigma);

public:
    EconomyParams ecoParams;
    vector < vector <double> > C;
    vector < vector <double> > eer;
    vector <double> Zexp;
    vector <double> normal_dist;
    int cur_t;
    int comp_ind;
};

bool pair_comp_fct(const pair<int,double> & a, const pair<int,double> & b)
{
    return a.second < b.second;
}

// Generate a random number following a normal distribution
double NationSave::getNormalDist(const double & mu, const double & sigma)
{
    const double u = (double) rand() / RAND_MAX;
    const double v = (double) rand() / RAND_MAX;
    const double r = sqrt( -2 * log(u) );
    const double x = r * sin(2 * PI * v);
    return sigma * x + mu;
}

double NationSave::find_zero1(EconomyParams & ecoParams, double WTm1,
                              double ZTm1, double x0, double x1, double tol)
{
    const double beta = ecoParams.beta;
    const double eta = ecoParams.eta;
    const double alpha = ecoParams.alpha;
    const double delta = ecoParams.delta;
    const double rho = ecoParams.rho;
    const double sigma = ecoParams.sigma;
    const int N = ecoParams.N;
    const int T = ecoParams.T;

    Zexp.resize(100);
    for (int i = 0; i < Zexp.size(); i++) {
        Zexp[i] = exp( rho * log(ZTm1) + getNormalDist(0.0, sigma) );
    }

    double err = 10.0;
    double x = (x0 + x1) / 2.0;
    while (err > tol) {
        x = (x0 + x1) / 2.0;

        double val = 0.0;
        for (int i = 0; i < Zexp.size(); i++) {
            const double v1 = pow( x / ( (1-delta) * (WTm1 - x) + Zexp[i] * pow(WTm1 - x, alpha) ), eta);
            val += beta * v1 * ( (1-delta) + Zexp[i] * alpha * pow(WTm1 - x, alpha-1) ) - 1;
        }
        val /= Zexp.size();

        err = abs(val);
        if (err <= tol) {
            break;
        } else if (val > 0) {
            x1 = x;
        } else {
            x0 = x;
        }
    }
    if (err > tol) {
        logout << "[WARNING] Could not reach precision: " << err << endl;
    }
    return x;
}


int NationSave::SetEconomyParameters(double beta, double eta, double alpha,
                        double delta, double rho, double sigma, int N, int T)
{
    logout << "\t\t" << "@ SetEconomyParameters: "
        << beta << ", " << eta << ", " << alpha << ", " << delta
        << ", " << rho << ", " << sigma << ", "
        << N << ", " << T
        << endl;

    ecoParams.beta = beta;
    ecoParams.eta = eta;
    ecoParams.alpha = alpha;
    ecoParams.delta = delta;
    ecoParams.rho = rho;
    ecoParams.sigma = sigma;
    ecoParams.N = N;
    ecoParams.T = T;

    srand(12345);

    C.assign(T, vector <double> (N, 0.0));
    eer.assign(T, vector <double> (N, 0.0));
    cur_t = -1;

    comp_ind = 0;

    return 0;
}

vector <double> NationSave::ConsumptionDecisionRule(vector <double> & Kt,
                                                    vector <double> & Zt)
{
    const double beta = ecoParams.beta;
    const double eta = ecoParams.eta;
    const double alpha = ecoParams.alpha;
    const double delta = ecoParams.delta;
    const double rho = ecoParams.rho;
    const double sigma = ecoParams.sigma;
    const int N = ecoParams.N;
    const int T = ecoParams.T;
    cur_t++;

    vector <double> Vt(N);
    vector <double> Wt(N);

    // Compute the wage
    for (int i = 0; i < N; i++) {
        Wt[i] = (1 - ecoParams.delta) * Kt[i] + Zt[i] * pow(Kt[i], ecoParams.alpha);
    }

    double kappa1 = 0.00001;
    if (sigma > 0.8 || rho > 0.75 || eta > 2.0) {
        kappa1 = 0.000001;
    }
#if USE_IMPROVED_VERSION
    double kappa2 = 0.5;
#else
    double kappa2 = 0.95;
#endif // USE_IMPROVED_VERSION

    // t == 0: Initial consumption
    if (cur_t == 0) {

        for (int i = 0; i < N; i++) {
            C[0][i] = kappa1 * Wt[i];
        }

    // t == T-1: Last time step
    } else if (cur_t == T-1) {

        // Simulate several values of Z_t in order to determine the consumption that
        //  minimizes the Euler residual
        for (int i = 0; i < N; i++) {
            Vt[i] = beta * (1 - delta + Zt[i] * alpha * pow(Kt[i], alpha - 1));
            double c1 = find_zero1(ecoParams, Wt[i], Zt[i], 0.0, Wt[i], 1e-10);
            C[cur_t][i] = c1;
        }

        // Compute the sum of the Euler residuals
        // Then, set the consumption of agent comp_ind to compensate for the non zero
        //  Euler residuals

#if USE_IMPROVED_VERSION
        double eer_tmp = 0.0;
        for (int i = 0; i < N; i++) {
            eer[cur_t][i] = pow( C[cur_t-1][i] / C[cur_t][i], eta ) * Vt[i] - 1.0;
            eer_tmp += eer[cur_t][i];
        }

        int trials_nb = 0;
        do {
            if (eer_tmp - eer[cur_t][comp_ind] < 0) {
                eer_tmp -= eer[cur_t][comp_ind];
                C[cur_t][comp_ind] = C[cur_t-1][comp_ind] * pow( Vt[comp_ind] / (1-eer_tmp), 1/eta );
                eer[cur_t][comp_ind] = pow( C[cur_t-1][comp_ind ] / C[cur_t][comp_ind], eta ) * Vt[comp_ind ] - 1.0;
                break;
            }
            trials_nb++;
            comp_ind++; // = rand() % N;
            if (comp_ind >= N) comp_ind = 0;
        } while (trials_nb < 1000);
#else
        double eer_tmp = 0.0;
        for (int i = 0; i < N; i++) {
            if (i != comp_ind) {
                eer[cur_t][i] = pow( C[cur_t-1][i] / C[cur_t][i], eta ) * Vt[i] - 1.0;
                eer_tmp += eer[cur_t][i];
            }
        }

        C[cur_t][comp_ind] = C[cur_t-1][comp_ind] * pow( Vt[comp_ind] / (1-eer_tmp), 1/eta );
        eer[cur_t][comp_ind] = pow( C[cur_t-1][comp_ind ] / C[cur_t][comp_ind ], eta ) * Vt[comp_ind ] - 1.0;
#endif // USE_IMPROVED_VERSION

    // Other time steps
    } else {

        for (int i = 0; i < N; i++) {
            Vt[i] = beta * (1 - delta + Zt[i] * alpha * pow(Kt[i], alpha - 1));
        }

        double eer_tmp = 0;
        for (int ii = 0; ii < N; ii++) {
            const int i = ii;

            // First possible consumption: a fraction of the wage
            double val_1 = kappa1 * Wt[i];
            // Second possible consumption: the optimal one
            double val_2 = pow(Vt[i], 1/eta) * C[cur_t-1][i];

            // Choose the highest consumption, unless it is higher than the wage
            if (val_1 >= val_2) {
                C[cur_t][i] = val_1;
            } else if (val_2 < kappa2 * Wt[i]) {
                C[cur_t][i] = val_2;
            } else {
                C[cur_t][i] = min( kappa2 * Kt[i], kappa2 * Wt[i] );
            }

            // Sanity check
            if (C[cur_t][i] == 0.0) {
                C[cur_t][i] = min( 0.05 * Kt[i], kappa2 * Wt[i] );
            }

            // Computation of the Euler residual
            eer[cur_t][i] = pow( C[cur_t-1][i] / C[cur_t][i], eta ) * Vt[i] - 1.0;
            eer_tmp += eer[cur_t][i];
        }

        // If the Euler residual os not zero, try to set it to zero
        if (abs(eer_tmp) > 1e-12) {
            if (eer_tmp > 0) {
            } else {
                // Try to find an agent which consumption can be set low enough for its
                //  Euler residual to compensate the other ones
                int trials_nb = 0;
                do {
                    comp_ind++;
                    if (comp_ind >= N) comp_ind = 0;
                    if (eer_tmp - eer[cur_t][comp_ind] < 0) {
                        eer_tmp -= eer[cur_t][comp_ind];
                        C[cur_t][comp_ind] = C[cur_t-1][comp_ind] * pow( Vt[comp_ind] / (1-eer_tmp), 1/eta );
                        eer[cur_t][comp_ind] = pow( C[cur_t-1][comp_ind ] / C[cur_t][comp_ind], eta ) * Vt[comp_ind ] - 1.0;
                        break;
                    }
                    trials_nb++;
                } while (trials_nb < 1000);
            }
        }
    }
    return C[cur_t];
}
