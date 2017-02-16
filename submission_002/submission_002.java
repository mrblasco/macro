import java.util.*;

public class NationSave {
    private double beta, eta, alpha, delta, factor, fpow, rho, fend;
    private int N, time, T;
    private double[] mult, pow;
    private List<Integer> order = new ArrayList<Integer>();
    private double[][] C, shocks, K, Z, W;
    private final Random rnd = new Random(99999);

    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
        this.beta = beta;
        this.eta = eta;
        this.alpha = alpha;
        this.delta = delta;
        this.rho = rho;
        this.N = N;
        this.T = T;
        C = new double[T + 1][N];
        Z = new double[T + 1][N];
        K = new double[T + 1][N];
        W = new double[T + 1][N];
        pow = new double[10000000];
        for (int i = 0; i < pow.length; i++) {
            pow[i] = Math.pow(i / 1000000.0, eta);
        }
        mult = new double[N];
        shocks = new double[T][N];
        for (int t = 0; t < T; t++) {
            for (int n = 0; n < N; n++) {
                shocks[t][n] = rnd.nextGaussian() * sigma;
            }
        }
        for (int n = 0; n < N; n++) {
            order.add(n);
        }
        factor = 1e-250;
        fpow = Math.exp(Math.log(1e-5 / factor) / T);
        return 0;
    }

    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt) {
        System.arraycopy(Kt, 0, K[time], 0, N);
        System.arraycopy(Zt, 0, Z[time], 0, N);
        final double[] Zt2 = Zt;
        Collections.sort(order, new Comparator<Integer>(){
            public int compare(Integer a, Integer b) {
                return Double.compare(Zt2[b], Zt2[a]);
            }
        });
        double fm = factor;
        if (time == T - 1) {
            double min = 9e99;
            double fmin = 0.0;
            double fmax = 0.7;
            double div = (fmax - fmin) / 100;
            for (int i = 0; i < 50; i++) {
                NEXT: for (double f = fmin; f <= fmax + div / 2 && f + div > f; f += div) {
                    double v = 0;
                    for (int j = 0; j < 60; j++) {
                        v += simulate(factor, f, 50, time, j);
                        if (v > min) continue NEXT;
                    }
                    if (v < min) {
                        min = v;
                        fend = f;
                    }
                }
                fmin = fend - div;
                fmax = fend + div;
                div /= 2;
            }
            fm += fend;
        }
        comp(Kt, Zt, fm, time, 1000);
        factor *= fpow;
        return C[time++];
    }

    private void comp(double[] Kt, double[] Zt, double fm, int t, int steps) {
        double[] Ct = C[t];
        double[] Wt = W[t];
        for (int n = 0; n < N; n++) {
            Wt[n] = (1 - delta) * Kt[n] + Zt[n] * Math.pow(Kt[n], alpha);
            Ct[n] = fm * Wt[n];
            mult[n] = beta * (1 - delta + Zt[n] * alpha * Math.pow(Kt[n], alpha - 1));
        }
        if (t > 0) {
            int cut = N / 5;
            double[] Ct1 = C[t - 1];
            double fmin = 1e-300;
            double fmax = 1;
            double best = 9e99;
            double max = t == T - 1 ? 1 : factor;
            double min = t == T - 1 ? 0 : -max;
            double div = (max - min) / 100;
            double fb = 0;
            for (int step = 0; step < steps; step++) {
                NEXT: for (double fi = min; fi <= max + div / 2 && fi + div > fi; fi += div) {
                    double err = 0;
                    for (int i = 0; i < order.size(); i++) {
                        int n = order.get(i);
                        double fn = fm + (i < cut ? -fi : fi);
                        if (fn <= fmin || fn >= fmax) continue NEXT;
                        err += pow(Ct1[n] / (fn * Wt[n])) * mult[n] - 1;
                    }
                    if (Math.abs(err) < best) {
                        best = Math.abs(err);
                        fb = fi;
                    }
                }
                min = fb - div;
                max = fb + div;
                if (max > 1) max = 1;
                if (min < 0 && t == T-1) min = 0;
                if (min < -1 && t < T-1) min = -1;
                div /= 2;
            }
            for (int i = 0; i < order.size(); i++) {
                int n = order.get(i);
                double fn = fm + (i < cut ? -fb : fb);
                Ct[n] = fn * Wt[n];
            }
        }
    }

    private double simulate(double f1, double f2, int steps, int t0, int q) {
        double ret = 0;
        for (int t = t0; t <= T; t++) {
            double[] st1 = t == 0 ? null : shocks[t - 1 - q];
            double[] Zt = Z[t];
            double[] Zt1 = t == 0 ? null : Z[t - 1];
            double[] Kt = K[t];
            double[] Wt = W[t];
            double[] Ct = C[t];
            double[] Kt1 = t == T ? null : K[t + 1];
            for (int n = 0; n < N; n++) {
                double Ktn = Kt[n];
                Wt[n] = (1 - delta) * Ktn + (t > t0 ? Zt[n] = Math.exp(rho * Math.log(Zt1[n]) + st1[n]) : Zt[n]) * Math.pow(Ktn, alpha);
            }
            if (t == T) System.arraycopy(Wt, 0, Ct, 0, N);
            else {
                if (f2 != 0 && t == T - 1) f1 += f2;
                if (f1 >= 1) return 1e99;
                comp(Kt, Zt, f1, t, steps);
                for (int n = 0; n < N; n++) {
                    Kt1[n] = Wt[n] - Ct[n];
                }
            }
            if (t > 0) {
                double e = 0;
                double[] Ct1 = C[t - 1];
                for (int n = 0; n < N; n++) {
                    e += beta * pow(Ct1[n] / Ct[n]) * (1 - delta + Zt[n] * alpha * Math.pow(Kt[n], alpha - 1)) - 1;
                }
                ret += Math.abs(e / N);
            }
        }
        return ret;
    }

    private double pow(double x) {
        double m = x * 1000000;
        int idx = (int) m;
        double f = m - idx;
        return idx >= pow.length - 1 ? Math.pow(x, eta) : pow[idx] * (1 - f) + pow[idx + 1] * f;
    }
}
