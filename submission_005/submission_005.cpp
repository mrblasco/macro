
#include <bits/stdc++.h>
using namespace std;

enum EconomyTypes {CNotInited, CIsNaN, CErrExceedsW, CErrTooSmall, CToSave, CToZero, CToRestore};

class Point
{
public:
    int idx;
    double Z;
    double C;

    Point () {}

    Point (int idx, double Z, double C)
    {
        this->idx = idx; this->Z = Z; this->C = C;
    }
};

class NationSave
{
public:
    const static bool printInfo = false;
    int test = 0;
    double beta, eta, alpha, delta, rho, sigma;
    int N, T;
    int currentT = 0;
    vector < vector <double> > W, C, K, Z; //first index - time, second index - subtest
    vector < vector <EconomyTypes> > E;

    double ScoreAnswer()
    {
        vector <double> errors(T);
        double eer = 0;
        for (int t = 0; t < T; t++)
        {
            double eer_t = 0;
            for (int i = 0; i < N; i++)
            {
                double eer_t_i = beta * pow(C[t][i] / C[t + 1][i], eta) *
                    (1 - delta + Z[t + 1][i] * alpha * pow(K[t + 1][i], alpha - 1)) - 1;
                eer_t += eer_t_i;
            }
            eer_t /= N;
            errors[t] = eer_t;
            if (t == (T - 1))
            {
                //Console.WriteLine(eer_t);
            }
            eer += abs(eer_t);
        }
        eer = 1 / (1 + eer / T);
        //if (this.test == 3)
        if (eer * 1000000 < 999500)
        {
            for (int i = 0; i < T; i++)
            {
                if (abs(errors[i]) > 1e-10)
                {
                    cout << "Time i " << i << " err " << errors[i] << endl;
                    for (int j = 0; j < N; j++)
                    {
                        double eer_t_i = beta * pow(C[i][j] / C[i + 1][j], eta) *
                             (1 - delta + Z[i + 1][j] * alpha *
                             pow(K[i + 1][j], alpha - 1)) - 1;

                        //if (Math.Abs(eer_t_i) > 0.5)
                        //    Console.WriteLine("Ct {0}  Ct+1 {1} err {2} n {3}",
                        //        C[i, j], C[i + 1, j], eer_t_i, j);
                    }
                }
            }
            cerr << eer << endl;
        }
        return eer * 1000000;
    }

    int SetEconomyParameters(double beta, double eta, double alpha, double delta,
        double rho, double sigma, int N, int T)
    {
        this->beta = beta; this->eta = eta; this->alpha = alpha;
        this->delta = delta; this->rho = rho; this->sigma = sigma;
        this->N = N; this->T = T;
        currentT = 0;
        W.assign(T+1, vector <double>(N));
        C.assign(T+1, vector <double>(N));
        K.assign(T+1, vector <double>(N));
        Z.assign(T+1, vector <double>(N));
        E.assign(T+1, vector <EconomyTypes> (N));
        for (int i = 0; i < N; i++)
        {
            K[0][i] = 1;
            Z[0][i] = 1;
        }
        if (printInfo)
        {
            cerr << "beta " << beta << " eta " << eta << " alpha " << alpha
                << " delta " << delta << " rho " << rho << " sigma " << sigma
                << " N " << N << " T " << T << endl;
        }
        return 1;
    }

    void AfterDecisionRule(vector <double> & Kt, vector <double> & Zt, vector <double> & Wt, vector <double> & Ct)
    {
        for (int i = 0; i < Zt.size(); i++)
        {
            K[currentT][i] = Kt[i];
            Z[currentT][i] = Zt[i];
            W[currentT][i] = Wt[i];
            C[currentT][i] = Ct[i];
            K[currentT + 1][i] = W[currentT][i] - C[currentT][i];
        }
        if (currentT == (T))
            FinEconomy();
        currentT++;
    }

    void FinEconomy()
    {
        for (int i = 0; i < N; i++)
        {
            C[T][i] = W[T][i];
        }
    }

    vector <double> ConsumptionDecisionRule(vector <double> Kt, vector <double> Zt)
    {
        int STEPS = 10000;
        double MAX_ERROR = 1e-300;
        double TO_RESTORE = 1e-250;
        double TO_RESTORE_MAX_COUNT = N/20;
        double AFTER_ERROR_K = 1e-3;
        double W_DIVIDE_TO = 10000;

        vector <double> Ct(Kt.size());
        vector <double> Wt(Kt.size());
        vector <double> Et(Kt.size());
        vector <Point> points(Kt.size());
        for (int i = 0; i < Ct.size(); i++)
        {
            Wt[i] = (1 - delta) * Kt[i] + Zt[i] * pow(Kt[i], alpha);
            points[i] = Point(i, Zt[i], (currentT > 0) ? C[currentT - 1][i] : 1);
            Ct[i] = Wt[i] / W_DIVIDE_TO;
        }
        if (true) {
            int cntCToZero = 0;
            int cntCToRestore = 0;
            vector <EconomyTypes> Ti(Kt.size());
            sort(points.begin(), points.end(),
                       [] (Point & r, Point & l) -> bool {
                            if (l.C < r.C) return false;
                            if (l.Z < r.Z) return false;
                            if (l.idx < r.idx) return false;
                            return true;
                       });
            int cnt = 0;
            for (int ii = 0; ii < points.size(); ii++)
            {
                int i = points[ii].idx;
                if (cnt >= Kt.size() / 2)
                {
                    Ti[i] = CToSave;
                }
                else
                {
                    if (currentT > 0 && currentT != (T-1) &&
                        C[currentT - 1][i] < TO_RESTORE && cntCToRestore < TO_RESTORE_MAX_COUNT)
                    {
                        Ti[i] = CToRestore;
                        cntCToRestore++;
                    }
                    else
                    {
                        Ti[i] = CToZero;
                        cntCToZero++;
                    }
                }
                cnt++;
            }
//*/
            for (int i = 0; i < Ct.size(); i++)
            {
                if (Ti[i] == CToSave)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                }
            }
            double CRateSum = 0;
            for (int i = 0; i < Ct.size(); i++)
            {
                if (Ti[i] == CToRestore && currentT >= 1)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                    Et[i] = beta * pow(C[currentT - 1][i] / Ct[i], eta) *
                                 (1 - delta + Zt[i] * alpha *
                                 pow(Kt[i], alpha - 1)) - 1;
                    {
                        CRateSum += Et[i];
                    }
                }
            }

            double sum_errors_plus = 0;
            int cnt_errors_plus = 0;
            double sum_errors_minus = Kt.size() / 2;
            for (int i = 0; i < Ct.size(); i++)
            {
                if (Ti[i] != CToSave)
                    continue;

                if (currentT == (T - 1))
                {
                    Ct[i] = Kt[i] / W_DIVIDE_TO;
                    double C1 = C[currentT - 1][i];
                    double C2_start = Wt[i] / 1.1;
                    double C2_end = C1 * 30;
                    double C2_best = -1;
                    double err_diff_best = std::numeric_limits<double>::max();
                    double err_best = std::numeric_limits<double>::max();

                    //double err_target = 1;
                    double err_target = (sum_errors_minus - sum_errors_plus) / (Kt.size() / 2 - cnt_errors_plus);

                    for (int s = 0; s <= STEPS; s++)
                    {

                        double C2 = C2_start + (C2_end - C2_start) * s / STEPS;
                        double Z3 = pow(Zt[i], rho);
                        double K3 = Wt[i] - C2;
                        double W3 = (1 - delta) * K3 + Z3 * pow(K3, alpha);
                        double C3 = W3;
                        double err = beta * pow(C2 / C3, eta) *
                                (1 - delta + Z3 * alpha *
                                pow(K3, alpha - 1)) - 1;
                        if (abs(err - err_target) < err_diff_best)
                        {
                            C2_best = C2;
                            err_diff_best = abs(err - err_target);
                            err_best = err;
                        }
                    }
                    if (C2_best > 0)
                    {
                        Ct[i] = C2_best;
                        Et[i] = beta * pow(C[currentT - 1][i] / Ct[i], eta) *
                             (1 - delta + Zt[i] * alpha *
                             pow(Kt[i], alpha - 1)) - 1;
                        sum_errors_plus += err_best;
                        cnt_errors_plus++;
                        CRateSum += Et[i];
                    }
                    continue;
                }

                if (currentT > 0)
                {
                    Et[i] = beta * pow(C[currentT - 1][i] / Ct[i], eta) *
                             (1 - delta + Zt[i] * alpha *
                             pow(Kt[i], alpha - 1)) - 1;
                    if (Et[i] < 0 && Et[i] > -1)
                    {
                        CRateSum += Et[i];
                    }
                    else
                    {
                        double curK = 1;
                        Ct[i] = C[currentT - 1][i] * pow(1 / curK, 1.0 / eta) *
                            pow(beta, 1.0 / eta) *
                            pow(1 - delta + Zt[i] * alpha * pow(Kt[i], alpha - 1), 1.0 / eta);
                        CRateSum += 0;
                    }
                    continue;
                }
            }

            for (int i = 0; i < Ct.size(); i++)
            {
                if (Ti[i] == CToZero)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                    if (currentT > 0)
                    {
                        double curK = 1 - CRateSum / cntCToZero;
                        Ct[i] = C[currentT - 1][i] * pow(1 / curK, 1.0 / eta) *
                            pow(beta, 1.0 / eta) *
                            pow(1 - delta + Zt[i] * alpha * pow(Kt[i], alpha - 1), 1.0 / eta);
                    }
                }
            }
        }
        {
            for (int i = 0; i < Ct.size(); i++)
            {
                if (std::isnan(Ct[i]))
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                }
                if (Ct[i] >= Wt[i] * (1 - MAX_ERROR))
                {
                    Ct[i] = Wt[i] * (1 - AFTER_ERROR_K);
                }
                if (Ct[i] < MAX_ERROR)
                {
                    Ct[i] = max(MAX_ERROR, Wt[i] * AFTER_ERROR_K);
                }
            }
            AfterDecisionRule(Kt, Zt, Wt, Ct);
            if (currentT == T - 1)
            {

            }
        }

        return Ct;
    }
};
