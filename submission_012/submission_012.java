import java.math.BigDecimal;
import java.util.*;

/**
 * Created by Nipuna on 8/24/2016.
 */
public class NationSave {
    static BigDecimal beta;
    static double betaDouble;
    static double deltaDouble;
    static double eta;
    static double alpha;
    static BigDecimal delta;
    static BigDecimal rho;
    static BigDecimal sigma;
    static Random rng = null;
    static int N;
    static int T;
    static double[] pC = null;
    static final BigDecimal startFraction = BigDecimal.valueOf(1.0).divide(BigDecimal.valueOf(10000.0));
    static int t;
    static double LIMIT = Math.pow(10, -200);
    static double LIMIT2 = Math.pow(10, -322);
    static int eater = -1;
    static double errSum = 0.0;
    static double[] errata = null;
    static int fixCount = 1;
    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
    {
        double[] C = new double[N];
        for (int i = 0; i < N; i++) {
            C[i] = getConsumption(Kt[i], Zt[i], i);
        }

        double eer_t = 0.0;
        if(t>0) {
            // eater
            boolean[] mark = new boolean[N];
            if (t == T-1)
            {
                Pair[] ps = new Pair[N];
                for (int i = 0; i < N; i++) {
                    ps[i] = new Pair(errata[i], i);
                }

                Arrays.sort(ps);
                eater = ps[N-1].index;
            }
            mark[eater] = true;
            for (int i = 0; i < N; i++) {
                if (mark[i])
                    continue;
                double eer_t_i = betaDouble * Math.pow(pC[i] / C[i], eta) * (1 - deltaDouble + Zt[i] * alpha * Math.pow(Kt[i], alpha - 1)) - 1;
                eer_t += eer_t_i;
            }
            C[eater] =    getEater(Kt[eater], Zt[eater], eater, eer_t);
        }

        pC = C;

        this.eater++;
        if (this.eater>=N)
            this.eater-=N;

        this.t++;
        return C;
    }


    private double getEater(double K, double Z, int simul, double amount) {
        BigDecimal Kb = BigDecimal.valueOf(K);
        BigDecimal Zb = BigDecimal.valueOf(Z);
        BigDecimal val = beta.multiply((BigDecimal.ONE.subtract(delta)).add((Zb.multiply(BigDecimal.valueOf(Math.pow(K, alpha-1)))).multiply(BigDecimal.valueOf(alpha))));
        BigDecimal W = ((BigDecimal.ONE.subtract(delta)).multiply(Kb)).add(Zb.multiply(BigDecimal.valueOf(Math.pow(K, alpha))));
        if (amount >= 1.0)
            amount = 0.0;
        double ret = (Math.pow(val.doubleValue()/(1.0 - amount), 1.0/eta) * pC[simul]);
        if (ret < LIMIT2 || W.doubleValue() < ret)
        {
            ret = startFraction.multiply(W).doubleValue();
        }
        return ret;

    }

    private double getConsumption(double K, double Z, int simul) {
        BigDecimal Kb = BigDecimal.valueOf(K);
        BigDecimal Zb = BigDecimal.valueOf(Z);
        BigDecimal W = ((BigDecimal.ONE.subtract(delta)).multiply(Kb)).add(Zb.multiply(BigDecimal.valueOf(Math.pow(K, alpha))));

        if(pC == null){
            return startFraction.multiply(W).doubleValue();
        }
        else{
            BigDecimal val = beta.multiply((BigDecimal.ONE.subtract(delta)).add((Zb.multiply(BigDecimal.valueOf(Math.pow(K, alpha-1)))).multiply(BigDecimal.valueOf(alpha))));
            double ret = (Math.pow(val.doubleValue(), 1.0/eta) * pC[simul]);

            if(t%200==199)
            {
                ret = startFraction.multiply(W).doubleValue();
            }
            if (t == T - 1)
            {
                ret = getR(K, Z, simul).multiply(W).doubleValue();
            }
            if(ret < LIMIT || ret > W.doubleValue())
            {
                ret = startFraction.multiply(W).doubleValue();
            }

            return ret;
        }
    }

    class Pair implements Comparable<Pair>{
        double value;
        int index;

        public Pair(double value, int index) {
            this.value = value;
            this.index = index;
        }

        @Override
        public int compareTo(Pair o) {
            if (this.value < o.value)
                return -1;
            else if (this.value > o.value)
                return 1;
            else
                return this.index-o.index;
        }
    }

    private BigDecimal getR(double K, double Z, int simul) {
        double w =  (1 - deltaDouble) * K + Z * Math.pow(K, alpha);;

        double bestR = 0.9;
        double r = 0.9;
        double minErr = Math.pow(10,200);
        int iter = 200;
        double step = (0.89)/(double)iter;
        int iterNormal = 1000;
        for (int i = 0; i < iter; i++) {
            double errTot = 0.0;
            for (int j = 0; j < iterNormal; j++) {
                double shock =  rng.nextGaussian() * sigma.doubleValue();
                double newZ = Math.exp(rho.doubleValue()*Math.log(Z) + shock);
                double newWage = (1 - deltaDouble) * w * (1.0 - r) + newZ * Math.pow(w * (1.0 - r), alpha);
                double err = (betaDouble*(Math.pow(w *r / newWage, eta) * (1.0 - deltaDouble + newZ * alpha * Math.pow(w * (1.0 - r), alpha - 1.0))) - (double)N/(double)(N-fixCount));
                errTot += err;
            }
            if (Math.abs(errTot) < minErr)
            {
                minErr = Math.abs(errTot);
                bestR = r;
            }
            r-=step;
        }
        BigDecimal ret = BigDecimal.valueOf(bestR);
        errata[simul] = Z;
        return ret;
    }

    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T){
        this.t = 0;
        this.pC = null;
        this.beta = BigDecimal.valueOf(beta);
        this.betaDouble = beta;
        this.eta = (eta);
        this.alpha = (alpha);
        this.delta = BigDecimal.valueOf(delta);
        this.deltaDouble = delta;
        this.rho = BigDecimal.valueOf(rho);
        this.sigma = BigDecimal.valueOf(sigma);
        this.N = N;
        this.T = T;
        this.eater = 0;
        this.errSum = 0.0;
        rng = new Random(1729);
        errata = new double[N];
        return 0;
    }
}
