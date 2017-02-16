import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class NationSave {
    private double beta;
    private double eta;
    private double alpha;
    private double delta;
    private double rho;
    private double sigma;
    private int N;
    private int T;

    private double[][] consumptionC;
    private double[] sumEER;

    private int t;

    private int count0;
    private int countGreaterThanWage;
    private int countNaN;

    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
        this.beta = beta;
        this.eta = eta;
        this.alpha = alpha;
        this.delta = delta;
        this.rho = rho;
        this.sigma = sigma;
        this.N = N;
        this.T = T;

        consumptionC = new double[T + 2][N];
        sumEER = new double[T + 2];

        t = 0;

        count0 = 0;
        countGreaterThanWage = 0;
        countNaN = 0;

        return 0;
    }

    public double[] ConsumptionDecisionRule(double[] capitalK, double[] shockZ) {
        ArrayList<Pair<Double, Integer>> list = new ArrayList<>();
        for (int n = 0; n < N; n++) {
            double wage = (1 - delta) * capitalK[n] + shockZ[n] * Math.pow(capitalK[n], alpha);
            list.add(new Pair<Double, Integer>(wage, n));
        }
        Collections.sort(list);

        for (int li = 0; li < list.size(); li++) {
            int n = list.get(li).second;

            double wage = (1 - delta) * capitalK[n] + shockZ[n] * Math.pow(capitalK[n], alpha);

            assert wage > 0 : Utils.toString("wage", wage, "(1 - delta) * capitalK[n]", (1 - delta) * capitalK[n], "shockZ[n] * Math.pow(capitalK[n], alpha)", shockZ[n] * Math.pow(capitalK[n], alpha), "capitalK[n]", capitalK[n]);

            double consumption;
            if (t == 0) {
                consumption = wage * 0.01;
            } else {
                if (li < 0.5 * N) {
                    consumption = wage * 0.01;
                } else {
                    double d = 1.0 - sumEER[t] / (N - li);
                    double e = beta * (1.0 - delta + shockZ[n] * alpha * Math.pow(capitalK[n], alpha - 1));
                    double x = Math.pow(d / e, 1.0 / eta);
                    consumption = consumptionC[t - 1][n] / x;
                }
            }

            if (consumption <= 1e-200) {
                consumption = wage * 0.01;
                count0++;
            }
            if (consumption >= wage) {
                consumption = wage * 0.01;
                countGreaterThanWage++;
            }
            if (Double.isNaN(consumption)) {
                consumption = wage * 0.01;
                countNaN++;
            }

            consumptionC[t][n] = consumption;

            if (t == 0) {
            } else {
                // double d = 1.0 - sumEER[t] / (N - li);

                double eer = eer(beta, eta, alpha, delta, consumptionC, capitalK, shockZ, t, n);
                sumEER[t] += eer;

                // Utils.debug("li", li, "sumEER[t]", sumEER[t], "eer", eer, "1.0 - sumEER[t] / (N - li)", d);
            }
        }

        if (t == T - 1) {
            Utils.debug("count0", count0, "countGreaterThanWage", countGreaterThanWage, "countNaN", countNaN);
        }

        double[] res = consumptionC[t];
        t++;

        return res;
    }

    private double eer(double beta, double eta, double alpha, double delta, double[][] C, double[] K, double[] Z, int t, int n) {
        return beta * Math.pow(C[t - 1][n] / C[t][n], eta) * (1.0 - delta + Z[n] * alpha * Math.pow(K[n], alpha - 1.0)) - 1.0;
    }

    public static void main(String[] args) {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(System.in))) {

            NationSave instance = new NationSave();

            int X = Integer.parseInt(reader.readLine());
            for (int x = 0; x < X; x++) {
                double beta = Double.parseDouble(reader.readLine());
                double eta = Double.parseDouble(reader.readLine());
                double alpha = Double.parseDouble(reader.readLine());
                double delta = Double.parseDouble(reader.readLine());
                double rho = Double.parseDouble(reader.readLine());
                double sigma = Double.parseDouble(reader.readLine());
                int N = Integer.parseInt(reader.readLine());
                int T = Integer.parseInt(reader.readLine());

                instance.SetEconomyParameters(beta, eta, alpha, delta, rho, sigma, N, T);
                for (int t = 0; t < T; t++) {
                    double[] Kt = new double[N];
                    for (int i = 0; i < N; i++) {
                        Kt[i] = Double.parseDouble(reader.readLine());
                    }
                    double[] Zt = new double[N];
                    for (int i = 0; i < N; i++) {
                        Zt[i] = Double.parseDouble(reader.readLine());
                    }

                    double[] Ct = instance.ConsumptionDecisionRule(Kt, Zt);
                    for (int i = 0; i < N; i++) {
                        System.out.println(Ct[i]);
                    }
                    System.out.flush();
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

class Pair<T extends Comparable<T>, S> implements Comparable<Pair<T, S>> {
    public T first;
    public S second;

    public Pair(T t, S s) {
        this.first = t;
        this.second = s;
    }

    private int hash = 0;

    @Override
    public int hashCode() {
        if (hash == 0) {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((first == null) ? 0 : first.hashCode());
            result = prime * result + ((second == null) ? 0 : second.hashCode());
            hash = result;
        }
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Pair<T, S> other = (Pair<T, S>) obj;
        if (first == null) {
            if (other.first != null)
                return false;
        } else if (!first.equals(other.first))
            return false;
        if (second == null) {
            if (other.second != null)
                return false;
        } else if (!second.equals(other.second))
            return false;
        return true;
    }

    @Override
    public int compareTo(Pair<T, S> o) {
        return first.compareTo(o.first);
    }
}

final class Utils {
    private Utils() {
    }

    public static final void debug(Object... o) {
        System.err.println(toString(o));
    }

    public static final void debug2(Object... o) {
        System.err.println(toString(o));
    }

    public static final String toString(Object... o) {
        return Arrays.deepToString(o);
    }

}
