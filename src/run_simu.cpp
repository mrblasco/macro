#include <bits/stdc++.h>
#include <sys/time.h>

// #include "NationSave.h"

using namespace std;

#define logout cerr
#define logouterr cerr
#define testerout cerr

class SimuEconomyParams
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

struct timeval tictoctime;
void tic()
{
    gettimeofday(&tictoctime, 0);
}
long int toc()
{
    struct timeval tictoctime2;
    gettimeofday(&tictoctime2, 0);
    return (tictoctime2.tv_sec - tictoctime.tv_sec) * 1e6 + (tictoctime2.tv_usec - tictoctime.tv_usec);
}
long int tocprint()
{
    struct timeval tictoctime2;
    gettimeofday(&tictoctime2, 0);
    const long int l = (tictoctime2.tv_sec - tictoctime.tv_sec) * 1e6 + (tictoctime2.tv_usec - tictoctime.tv_usec);
    std::logout << "Elapsed time: " << (float)l / 1e6 << std::endl;
    return l;
}

int main(int argc, char* argv[])
{

    logout << "*** Files ***" << endl;

    string params_filename = "./economies.csv";
    if (argc > 1) params_filename = argv[1];
    cout << "\t" << "params_filename: " << params_filename << endl;
    string shocks_filename = "./shocks_system.csv";
    if (argc > 2) shocks_filename = argv[2];
    cout << "\t" << "shocks_filename: " << shocks_filename << endl;
    string output_filename = "./albantor30.csv";
    if (argc > 3) output_filename = argv[3];
    cout << "\t" << "output_filename: " << output_filename << endl;
    string log_filename = output_filename + ".log";
    int test_offset = 0;
    if (argc > 4) test_offset = atoi(argv[4]);
    cout << "\t" << "test_offset: " << test_offset << endl;
    int num_tests = -1;
    if (argc > 5) num_tests = atoi(argv[5]);
    cout << "\t" << "num_tests: " << num_tests << endl;

    logout << "*** Loading simulation parameters ***" << endl;
    tic();

    const int N = 50;
    const int T = 2000;
    const double K0 = 1.0;
    const double Z0 = 1.0;

    vector <SimuEconomyParams> ecoParams;
    {
        ifstream fin(params_filename);
        if (! fin.good()) {
            cerr << "[ERROR] Could not open economies file: " << params_filename << endl;
            return -1;
        }
        string line;
        getline(fin, line);
        while(getline(fin, line)) {
            if (line.back() == '\r') line.pop_back();
            if (line.empty()) continue;

            stringstream ss(line);
            char c;
            SimuEconomyParams params;
            ss >> params.beta >> c
                >> params.eta >> c
                >> params.alpha >> c
                >> params.delta >> c
                >> params.rho >> c
                >> params.sigma;
            params.N = N;
            params.T = T;
            ecoParams.push_back(params);
        }
        fin.close();
    }
    int numTestCases = ecoParams.size();
    logout << "\t" << "numTestCases: " << numTestCases << endl;

    cout << "\t" << "Settings: "
        << "N: " << N << ", T: " << T
        << ", K0: " << K0 << ", Z0: " << Z0 << endl;

    logout << "\t"; tocprint();


    logout << "*** Loading shocks ***" << endl;
    tic();

    vector < vector <double> > shocks(N, vector <double> (T));
    {
        ifstream fin(shocks_filename);
        if (! fin.good()) {
            cerr << "[ERROR] Could not open shocks file: " << shocks_filename << endl;
            return -1;
        }
        for (int t = 0; t < T; t++) {
            for (int n = 0; n < N; n++) {
                fin >> shocks[n][t];
                if (n < N-1) fin.ignore(1);
            }
        }
        fin.close();
    }

    logout << "\t"; tocprint();


    logout << "*** Running simulation ***" << endl;
    tic();

    // Output file
    ofstream fout(output_filename);
    if (! fout.good()) {
        cerr << "[ERROR] Could not open output file: " << output_filename << endl;
        return -1;
    }
    fout << setprecision(17);
    fout << "# Economy_id, simulation_id, period, C, Z, K, EulerResidual" << endl;

    // Log file
    ofstream fout_log(log_filename);
    if (! fout_log.good()) {
        cerr << "[ERROR] Could not open log file: " << log_filename << endl;
        return -1;
    }
    fout_log << setprecision(10);
    fout_log << "# Economy_id, score" << endl;

    double totalScore = 0.0;
    int last_test = numTestCases;
    if (num_tests >= 0) last_test = min(numTestCases, test_offset + num_tests);
    for (int test = test_offset; test < last_test; test++) {
        logout << "\t" << "case: " << setw(6) << test+1 << " / " << last_test << endl;

        NationSave predictor;

        vector < vector <double> > C(T + 1, vector <double> (N, 0.0));
        vector < vector <double> > K(T + 1, vector <double> (N, K0));
        vector < vector <double> > W(T + 1, vector <double> (N, 0.0));
        vector < vector <double> > Z(T + 1, vector <double> (N, Z0));
        vector < vector <double> > eer(T + 1, vector <double> (N, 0.0));

        bool hasError = false;

        predictor.SetEconomyParameters(ecoParams[test].beta, ecoParams[test].eta, ecoParams[test].alpha, ecoParams[test].delta, ecoParams[test].rho, ecoParams[test].sigma, ecoParams[test].N, ecoParams[test].T);

        for (int t = 1; t <= ecoParams[test].T; t++) {
            for (int n = 0; n < ecoParams[test].N; n++) {
                Z[t][n] = exp(ecoParams[test].rho * log( Z[t - 1][n] ) + ecoParams[test].sigma * shocks[n][t-1]);
            }
        }

        for (int t = 0; t <= ecoParams[test].T; t++) {

            for (int n = 0; n < ecoParams[test].N; n++) {
                // if (t > 0)  Z[t][n] = exp(ecoParams[test].rho * log( Z[t - 1][n] ) + shocks[n][t-1]);
                W[t][n] = (1 - ecoParams[test].delta) * K[t][n] + Z[t][n] * pow(K[t][n], ecoParams[test].alpha);
            }

            if (t != ecoParams[test].T) {
                C[t] = predictor.ConsumptionDecisionRule(K[t], Z[t]);

                for (int n = 0; n < ecoParams[test].N; n++) {
                    // Sanity check
                    if (C[t][n] < 0.0) {
                        logouterr << "[ERROR] C[" << t << "][" << n << "] cannot be negative." << endl;
                        ecoParams[test].score = -1;
                        hasError = true;
                        break;
                    }
                    if (C[t][n] > W[t][n]) {
                        logouterr << "[ERROR] C[" << t << "][" << n << "] cannot exceed wage." << endl;
                        ecoParams[test].score = -2;
                        hasError = true;
                        break;
                    }
                    if (std::isnan(C[t][n])) {
                        logouterr << "[ERROR] C[" << t << "][" << n << "] must be a real number." << endl;
                        logouterr << "\t" << C[0][n] << endl;
                        ecoParams[test].score = -3;
                        hasError = true;
                        break;
                    }
                    // Update capital for next period
                    K[t + 1][n] = W[t][n] - C[t][n];
                }

            } else {
                // On last time period, consume all of the wage
                C[t] = W[t];
            }
/*
            if (t < 3) {
                if (t == 0) {
                    cout << "[" << setw(4) << "t" << "]"
                        << " \t" << setw(20) << "K[t][n]"
                        << " \t" << setw(20) << "Z[t][n]"
                        << " \t" << setw(20) << "W[t][n]"
                        << " \t" << setw(20) << "C[t][n]"
                        << endl;
                }
                for (int n = 0; n < 4; n++) {
                    cout << setprecision(16);
                    cout << "[" << setw(4) << t << "]"
                        << " \t" << setw(20) << K[t][n]
                        << " \t" << setw(20) << Z[t][n]
                        << " \t" << setw(20) << W[t][n]
                        << " \t" << setw(20) << C[t][n]
                        << endl;
                }
            }
//*/
            if (hasError) break;

            if (t > 0) {
                for (int n = 0; n < ecoParams[test].N; n++) {
                    eer[t-1][n] = ecoParams[test].beta * pow(C[t-1][n] / C[t][n], ecoParams[test].eta) * (1 - ecoParams[test].delta + Z[t][n] * ecoParams[test].alpha * pow(K[t][n], ecoParams[test].alpha - 1)) - 1;
                }
            }

            if (hasError) break;
        }

        if (! hasError) {
            ecoParams[test].score = 0;
            for (int t = 0; t < ecoParams[test].T; t++) {
                double eer_t = 0;
                for (int n = 0; n < ecoParams[test].N; n++) {
                    eer_t += eer[t][n];
                }
                ecoParams[test].score += abs(eer_t) / ecoParams[test].N;
            }
            ecoParams[test].score = 1e6 / (1 + ecoParams[test].score / ecoParams[test].T);
            totalScore += ecoParams[test].score;
        }

        cout << "\t\t" << "score: " << setprecision(10) << ecoParams[test].score << endl;

        // Output simulation data
        if (! (ecoParams[test].score < 0)) {
            cout << "\t\t" << "Saving..." << endl;
            // fout << "# Economy_id, simulation_id, period, C, Z, K, EulerResidual" << endl;
            for (int n = 0; n < N; n++) {
                for (int t = 0; t <= T; t++) {
                    fout << (test) << "," << n << "," << t
                        << "," << C[t][n] << "," << Z[t][n] << "," << K[t][n] << ",";
                    if (t < T) fout << eer[t][n];
                    fout << endl;
                }
            }
            if (! fout.good()) {
                cerr << "[ERROR] Error writing to output file." << endl;
                return -1;
            }
        }
        // Output log info
        fout_log << test << ",\t" << ecoParams[test].score << endl;

    }
    totalScore /= numTestCases;

    fout.close();
    fout_log.close();

    cout << endl
        << "\t" << "Total score: " << totalScore << endl;


    logout << "\t"; tocprint();

    return 0;
}
