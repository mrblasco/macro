using System;
using System.Collections.Generic;

public enum EconomyTypes {CNotInited, CIsNaN, CErrExceedsW, CErrTooSmall, CToSave, CToZero, CToRestore};

public class NationSave
{
    static bool printInfo = false;
    public int test = 0;
    double beta, eta, alpha, delta, rho, sigma;
    int N, T;
    int currentT = 0;
    public double[,] W, C, K, Z; //first index - time, second index - subtest
    public EconomyTypes[,] E;

    public double ScoreAnswer()
    {
        double[] errors = new double[T];
        double eer = 0;
        for (int t = 0; t < T; t++)
        {
            double eer_t = 0;
            for (int i = 0; i < N; i++)
            {
                double eer_t_i = beta * Math.Pow(C[t, i] / C[t + 1, i], eta) * 
                    (1 - delta + Z[t + 1, i] * alpha * Math.Pow(K[t + 1, i], alpha - 1)) - 1;
                eer_t += eer_t_i;
            }
            eer_t /= N;
            errors[t] = eer_t;
            if (t == (T - 1))
            {
                //Console.WriteLine(eer_t);
            }
            eer += Math.Abs(eer_t);
        }
        eer = 1 / (1 + eer / T);
        //if (this.test == 3)
        if (eer * 1000000 < 999500)
        {
            for (int i = 0; i < T; i++)
            {
                if (Math.Abs(errors[i]) > 1e-10)
                { 
                    Console.WriteLine("Time i {0} err {1}", i, errors[i]);
                    for (int j = 0; j < N; j++)
                    {
                        double eer_t_i = beta * Math.Pow(C[i, j] / C[i + 1, j], eta) *
                             (1 - delta + Z[i + 1, j] * alpha * 
                             Math.Pow(K[i + 1, j], alpha - 1)) - 1;

                        //if (Math.Abs(eer_t_i) > 0.5)
                        //    Console.WriteLine("Ct {0}  Ct+1 {1} err {2} n {3}", 
                        //        C[i, j], C[i + 1, j], eer_t_i, j);
                    }
                }
            }
            Console.Error.WriteLine(eer);
        }
        return eer * 1000000;
    }
  
    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, 
        double rho, double sigma, int N, int T)
    {
        this.beta = beta; this.eta = eta; this.alpha = alpha;
        this.delta = delta; this.rho = rho; this.sigma = sigma;
        this.N = N; this.T = T;
        currentT = 0;
        W = new double[T + 1, N];
        C = new double[T + 1, N];
        K = new double[T + 1, N];
        Z = new double[T + 1, N];
        E = new EconomyTypes[T + 1, N];
        for (int i = 0; i < N; i++)
        {
            K[0, i] = 1;
            Z[0, i] = 1;
        }
        if (printInfo)
        { 
            Console.Error.WriteLine(
                "beta {0} eta {1} alpha {2} delta {3} rho {4} sigma {5} N {6} T {7}",
                beta, eta, alpha, delta, rho, sigma, N, T);
        }
        return 1;
    }
    
    public void AfterDecisionRule(double[] Kt, double[] Zt, double[] Wt, double[] Ct)
    {
        for (int i = 0; i < Zt.Length; i++)
        {            
            K[currentT, i] = Kt[i];
            Z[currentT, i] = Zt[i];
            W[currentT, i] = Wt[i];
            C[currentT, i] = Ct[i];
            K[currentT + 1, i] = W[currentT, i] - C[currentT, i];
        }
        if (currentT == (T))
            FinEconomy();
        currentT++;
    }

    public void FinEconomy()
    {
        for (int i = 0; i < N; i++)
        {
            C[T, i] = W[T, i];
        }
    }

    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
    {
        int STEPS = 10000;
        double MAX_ERROR = 1e-300;
        double TO_RESTORE = 1e-250;
        double TO_RESTORE_MAX_COUNT = N/20;
        double AFTER_ERROR_K = 1e-3;
        double W_DIVIDE_TO = 10000;

        double[] Ct = new double[Kt.Length];
        double[] Wt = new double[Kt.Length];
        double[] Et = new double[Kt.Length];
        Point[] points = new Point[Kt.Length];
        for (int i = 0; i < Ct.Length; i++)
        {
            Wt[i] = (1 - delta) * Kt[i] + Zt[i] * Math.Pow(Kt[i], alpha);
            points[i] = new Point(i, Zt[i], (currentT > 0) ? C[currentT - 1, i] : 1);
            Ct[i] = Wt[i] / W_DIVIDE_TO;             
        }
        try
        {
            int cntCToZero = 0;
            int cntCToRestore = 0;
            EconomyTypes[] Ti = new EconomyTypes[Kt.Length];
            Array.Sort(points);
            int cnt = 0;
            foreach (Point p in points)
            {
                int i = p.idx;
                if (cnt >= Kt.Length / 2)
                {
                    Ti[i] = EconomyTypes.CToSave;
                }
                else
                {
                    if (currentT > 0 && currentT != (T-1) &&
                        C[currentT - 1, i] < TO_RESTORE && cntCToRestore < TO_RESTORE_MAX_COUNT)
                    {
                        Ti[i] = EconomyTypes.CToRestore;
                        cntCToRestore++;
                    }
                    else
                    {
                        Ti[i] = EconomyTypes.CToZero;
                        cntCToZero++;
                    }
                }
                cnt++;
            }
            for (int i = 0; i < Ct.Length; i++)
            {
                if (Ti[i] == EconomyTypes.CToSave)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                }
            }
            double CRateSum = 0;
            for (int i = 0; i < Ct.Length; i++)
            {
                if (Ti[i] == EconomyTypes.CToRestore && currentT >= 1)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                    Et[i] = beta * Math.Pow(C[currentT - 1, i] / Ct[i], eta) *
                                 (1 - delta + Zt[i] * alpha *
                                 Math.Pow(Kt[i], alpha - 1)) - 1;
                    {
                        CRateSum += Et[i];
                    }
                }
            }

            double sum_errors_plus = 0;
            int cnt_errors_plus = 0;
            double sum_errors_minus = Kt.Length / 2;
            for (int i = 0; i < Ct.Length; i++)
            {
                if (Ti[i] != EconomyTypes.CToSave)
                    continue;

                if (currentT == (T - 1))
                {
                    Ct[i] = Kt[i] / W_DIVIDE_TO;
                    double C1 = C[currentT - 1, i];
                    double C2_start = Wt[i] / 1.1;
                    double C2_end = C1 * 30;
                    double C2_best = -1;
                    double err_diff_best = double.MaxValue;
                    double err_best = double.MaxValue;

                    //double err_target = 1;
                    double err_target = (sum_errors_minus - sum_errors_plus) / (Kt.Length / 2 - cnt_errors_plus);

                    for (int s = 0; s <= STEPS; s++)
                    {
                        
                        double C2 = C2_start + (C2_end - C2_start) * s / STEPS;
                        double Z3 = Math.Pow(Zt[i], rho);
                        double K3 = Wt[i] - C2;
                        double W3 = (1 - delta) * K3 + Z3 * Math.Pow(K3, alpha);
                        double C3 = W3;
                        double err = beta * Math.Pow(C2 / C3, eta) *
                                (1 - delta + Z3 * alpha *
                                Math.Pow(K3, alpha - 1)) - 1;
                        if (Math.Abs(err - err_target) < err_diff_best)
                        {
                            C2_best = C2;
                            err_diff_best = Math.Abs(err - err_target);
                            err_best = err;
                        }
                    }
                    if (C2_best > 0)
                    {
                        Ct[i] = C2_best;
                        Et[i] = beta * Math.Pow(C[currentT - 1, i] / Ct[i], eta) *
                             (1 - delta + Zt[i] * alpha *
                             Math.Pow(Kt[i], alpha - 1)) - 1;
                        sum_errors_plus += err_best;
                        cnt_errors_plus++;
                        CRateSum += Et[i];
                    }
                    continue;
                }

                if (currentT > 0)
                {
                    Et[i] = beta * Math.Pow(C[currentT - 1, i] / Ct[i], eta) *
                             (1 - delta + Zt[i] * alpha *
                             Math.Pow(Kt[i], alpha - 1)) - 1;
                    if (Et[i] < 0 && Et[i] > -1)
                    {
                        CRateSum += Et[i];
                    }
                    else
                    {
                        double curK = 1;
                        Ct[i] = C[currentT - 1, i] * Math.Pow(1 / curK, 1.0 / eta) *
                            Math.Pow(beta, 1.0 / eta) *
                            Math.Pow(1 - delta + Zt[i] * alpha * Math.Pow(Kt[i], alpha - 1), 1.0 / eta);
                        CRateSum += 0;
                    }
                    continue;
                }
            }

            for (int i = 0; i < Ct.Length; i++)
            {
                if (Ti[i] == EconomyTypes.CToZero)
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                    if (currentT > 0)
                    {
                        double curK = 1 - CRateSum / cntCToZero;
                        Ct[i] = C[currentT - 1, i] * Math.Pow(1 / curK, 1.0 / eta) *
                            Math.Pow(beta, 1.0 / eta) *
                            Math.Pow(1 - delta + Zt[i] * alpha * Math.Pow(Kt[i], alpha - 1), 1.0 / eta);
                    }
                }
            }
        }
        catch (Exception e)
        {
            Console.Error.WriteLine(e);
        }
        finally
        {
            for (int i = 0; i < Ct.Length; i++)
            {
                if (Double.IsNaN(Ct[i]))
                {
                    Ct[i] = Wt[i] / W_DIVIDE_TO;
                }
                if (Ct[i] >= Wt[i] * (1 - MAX_ERROR))
                {
                    Ct[i] = Wt[i] * (1 - AFTER_ERROR_K);
                }
                if (Ct[i] < MAX_ERROR)
                {
                    Ct[i] = Math.Max(MAX_ERROR, Wt[i] * AFTER_ERROR_K);
                }
            }
            AfterDecisionRule(Kt, Zt, Wt, Ct);
            if (currentT == T - 1)
            {

            }
        }
        return Ct;
    }
}


public class Point : IComparable<Point>
{
    public int idx;
    public double Z;
    public double C;
    
    public Point (int idx, double Z, double C)
    {
        this.idx = idx; this.Z = Z; this.C = C;
    }

    int IComparable<Point>.CompareTo(Point other)
    {
        if (other.C > this.C) return 1;
        if (other.Z > this.Z) return 1;
        if (other.idx > this.idx) return 1;
        return 0;
    }
}
