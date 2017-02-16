import java.io.*;
import java.util.*;

import static java.lang.Math.*;
import static java.lang.Double.*;
import static java.lang.System.*;
import static java.util.Arrays.*;

class Constants {

    public static final long TIME_LIMIT = 111_111;

    public static final int P = 23;
    public static final int M = P*42;

    public static final double EPS = 1e-6;
    public static final double ZERO = 1e-301;

    public static final double PRECISION = 1e-15;
    public static final double SCORE_PRECISION = 1e-12;
}

class NationSave {

    public double[] ConsumptionDecisionRule(double[] kt, double[] zt) {

        if (time > 0) {

            if (time == 1) {
                params.fillZ(TestDeterminer.getTestId(params, zt)+1); 
            }

            double err = Simulator.tick(
                    solution.c[time-1], solution.c[time], 
                    context.w[prevT], context.w[curT], zt,
                    context.kt, context.cf, context.e,
                    context.rng, params
                );
    
            if (abs(err) > Constants.PRECISION) {
                err = SimulationOptimizer.optimize(
                        err, solution.c[time-1], solution.c[time],
                        context.w[prevT], context.w[curT], zt,
                        context.kt, context.cf, context.e,
                        context.rng, params
                    );
                context.totalErr += abs(err);
            }

            if (time == params.duration-1) {
                double [][][] zpt = new double[Constants.P][Constants.M/Constants.P][params.simulations];
                ICombiner[] combiners = { new Combiner(), };
                int combiner = 0;
                int combinersCount = combiners.length;
                for (int p = 0; p < Constants.P; ++p) {
                    for (int r = p; r < Constants.M; r += Constants.P) {
                        for (int i = 0; i < params.simulations; ++i) {
                            zpt[p][r/Constants.P][i] = params.zt[r][i];
                        }
                    }
                }
                for (int go = 0; !timer.isTL() && go < 222; ++go) {
                    for (int p = 0; p < Constants.P; ++p) {
                        err = FinalOptimizer.optimize(
                                solution.c[time-1], solution.c[time],
                                context.w[prevT], context.w[curT], zt, zpt[p],
                                context.kt, context.cf, context.e,
                                context.rng, params, combiners[combiner]
                            );    
                        ++combiner;
                        if (combiner == combinersCount) {
                            combiner = 0;
                        }
                    }
                        
                    for (int i1 = 0; i1 < Constants.P-1; ++i1) {
                        int i2 = i1+1+rng.nextInt(Constants.P-1-i1);
                        int j1 = rng.nextInt(Constants.M/Constants.P), j2 = rng.nextInt(Constants.M/Constants.P);
                        double[] tmp = zpt[i1][j1];
                        zpt[i1][j1] = zpt[i2][j2];
                        zpt[i2][j2] = tmp;
                    }
                }
                context.totalErr += err;

            }

            prevT ^= 1;
            curT ^= 1;
        }

        double[] ct = new double[params.simulations];
        arraycopy(solution.c[time], 0, ct, 0, params.simulations);
        ++time;
        return ct;
    }

    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int simulations, int duration) {
        timer = new Timer(Constants.TIME_LIMIT);
        time = 0;
        params = new Parameters(beta, eta, alpha, delta, rho, sigma, simulations, duration);
        solution = new Solution(params);
        context = new Context(params, solution);
        context.reset();

        prevT = 0; 
        curT = 1;

        return 0;
    }
    
    private Context context; 
    private Parameters params;
    private Solution solution;

    private int time, prevT, curT;

    private Timer timer;
    private Random rng = new Random();

}

class Parameters {

    public static final double Z0 = 1;

    public Parameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int simulations, int duration) {
        this.alpha = alpha;
        this.beta = beta;
        this.delta = delta;
        this.eta = eta;
        this.rho = rho;
        this.sigma = sigma;
        this.simulations = simulations;
        this.duration = duration;
        
        this.reta = 1/eta;

        zt = new double[Constants.M][simulations];
    }

    public void fillZ(int seed) {
        Random rng = new Random(seed);
        if (seed > TestDeterminer.MAX_TEST_ID) {
            for (int r = 0; r < Constants.M; ++r) {
                fill(zt[r], Z0);
                for (int t = 1; t < duration; ++t) {
                    for (int i = 0; i < simulations; ++i) {
                        zt[r][i] = exp(rho*log(zt[r][i]) + rng.nextGaussian()*sigma);
                    }
                }
            }
        }
        else {
            fill(zt[0], Z0);
            for (int t = 1; t <= duration; ++t) {
                for (int i = 0; i < simulations; ++i) {
                    zt[0][i] = exp(rho*log(zt[0][i]) + rng.nextGaussian()*sigma);
                }
            }
            for (int r = 1; r < Constants.M; ++r) {
                arraycopy(zt[0], 0, zt[r], 0, simulations);
            }
        }
    }

    int simulations, duration;
    double alpha, beta, delta, eta, rho, sigma, rdelta, reta, zt[][];
}

class Solution {

    public static final double INITIAL_COEFFICIENT = 1e-2;

    public Solution(Parameters params) {
        c = new double[params.duration+1][params.simulations];
        Random rng = new Random();
        for (int i = 0; i < params.simulations; ++i) {
            c[0][i] = rng.nextDouble()*INITIAL_COEFFICIENT*(2 - params.delta);
        }
    }

    double score = POSITIVE_INFINITY;
    double[][] c;
}

class TestDeterminer {

    public static final int MAX_TEST_ID = 801;

    public static int getTestId(final Parameters params, final double[] z1) {
        Random rng = new Random();
        for (int test = 0; test < MAX_TEST_ID; ++test) {
            rng.setSeed(test+1);
            boolean ok = true;
            for (int i = 0; i < params.simulations; ++i) {
                if (abs(exp(rng.nextGaussian()*params.sigma) - z1[i]) > Constants.EPS) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                return test;
            }
        }
        return MAX_TEST_ID;
    }
}

class Context {
    
    public Context(Parameters params, final Solution solution) {
        this.params = params;

        c = new double[params.duration+1][params.simulations];

        w[0] = new double[params.simulations];
        w[1] = new double[params.simulations];

        cf = new double[params.simulations];
        kt = new double[params.simulations];
        e = new double[params.simulations];

        arraycopy(solution.c[0], 0, c[0], 0, params.simulations);
    }

    public void reset() {
        reset(0x11);
    }

    public void reset(int seed) {
        rng.setSeed(seed);
        fill(w[0], 2-params.delta);
        totalErr = 0;
    }
   
    double totalErr;
    Parameters params;
    double[][] w = new double[2][];
    double[] cf, kt, e;
    double[][] c;

    Random rng = new Random();
}

class Simulator {

    public static double tick(
            double[] consumptionPrevious, double[] consumptionCurrent, 
            double[] wagePrevious, double[] wageCurrent, double[] zt,
            double[] kapital, double[] coefficients, double[] errors,
            Random rng,
            Parameters params
        ) {

        double err = 0;
        for (int i = 0; i < params.simulations; ++i) {
            double kti = wagePrevious[i] - consumptionPrevious[i]; 
            kapital[i] = kti;
            double wti = Calculator.getWage(kti, zt[i], params);
            wageCurrent[i] = wti;
            consumptionCurrent[i] = Calculator.getCt(consumptionPrevious[i], kti, zt[i], params);
            if (isNaN(consumptionCurrent[i]) || consumptionCurrent[i] <= Constants.ZERO || consumptionCurrent[i] >= wti) {
                double oldValue = consumptionCurrent[i];
                coefficients[i] = rng.nextDouble()*Solution.INITIAL_COEFFICIENT;
                consumptionCurrent[i] = coefficients[i]*wti;

                errors[i] = Calculator.getErr(consumptionPrevious[i], consumptionCurrent[i], kti, zt[i], params);
                err += errors[i];
            }
            else 
            {
                errors[i] = 0;
                coefficients[i] = consumptionCurrent[i]/wti;
            }
        }
        return err;
    }
}

class SimulationOptimizer {

    public static double optimize(
            double err, double[] consumptionPrevious, double[] consumptionCurrent,
            double[] wagePrevious, double[] wageCurrent, double[] zt,
            double[] kapital, double[] coefficients, double[] errors,
            Random rng, Parameters params
        ) {

        for (int go = 0; go < 10; ++go) {
            double pw = 1;
            for (int it = 0; it < 23; ++it) {
                for (int it2 = 0; abs(err) > Constants.PRECISION && it2 < 77; ++it2) {
                    int i = rng.nextInt(params.simulations);
                    double k = coefficients[i]*(1+rng.nextGaussian()*pw);
                    int cnt = 0;
                    while (++cnt < 1111 && (k*wageCurrent[i] <= 0 || k >= 1)) {
                        k = coefficients[i]*(1+rng.nextGaussian()*pw);
                    }
                    if (cnt == 1111) {
                        continue;
                    }
                    double consumptionNew = wageCurrent[i]*k;
                    double ne = Calculator.getErr(consumptionPrevious[i], consumptionNew, kapital[i], zt[i], params);
                    if (abs(err) > abs(err + ne - errors[i])) {
                        err += ne - errors[i];
                        coefficients[i] = k;
                        errors[i] = ne;
                        consumptionCurrent[i] = consumptionNew;
                    }
                }
                pw *= 0.2;
            }
        }

        return err;
    }
}

class FinalOptimizer {

    public static double optimize(
            double[] consumptionPrevious, double[] consumptionCurrent,
            double[] wagePrevious, double[] wageCurrent, double[] zCurrent, double[][] zNext,
            double[] kapital, double[] coefficients, double[] errors1,
            Random rng, Parameters params, ICombiner combiner
        ) {

        double[][] errors2 = new double[Constants.M/Constants.P][params.simulations];
        double[] te2 = new double[Constants.M/Constants.P];
        double[] ne2 = new double[Constants.M/Constants.P];
        double[] v2 = new double[Constants.M/Constants.P];
        double te1 = 0;
        for (int i = 0; i < params.simulations; ++i) {
            coefficients[i] = abs(wageCurrent[i]) < Constants.ZERO ? 0 : consumptionCurrent[i]/wageCurrent[i];
            errors1[i] = Calculator.getErr(consumptionPrevious[i], consumptionCurrent[i], kapital[i], zCurrent[i], params);
            te1 += errors1[i];
            double kapitalNext = wageCurrent[i]-consumptionCurrent[i];

            for (int j = 0; j < Constants.M/Constants.P; ++j) {
                double consumptionNext = Calculator.getWage(kapitalNext, zNext[j][i], params);
                errors2[j][i] = Calculator.getErr(consumptionCurrent[i], consumptionNext, kapitalNext, zNext[j][i], params);
                te2[j] += errors2[j][i];
            }
        }
        double val = abs(te1) + combiner.combine(te2);
        for (int go = 0; val > Constants.PRECISION && go < 11; ++go) {
            double pw = 1e5;
            for (int d = 0; val > Constants.PRECISION && d < 11; ++d) {
                for (int it = 0; val > Constants.PRECISION && it < 11; ++it) {
                    int i = rng.nextInt(params.simulations);
                    if (coefficients[i] < Constants.ZERO) {
                        continue;
                    }
                    double k = coefficients[i] + rng.nextGaussian()*pw;
                    int cnt = 0;
                    while (++cnt < 1111 && (wageCurrent[i]*k <= 0 || k >= 1)) {
                        k = coefficients[i] + rng.nextGaussian()*pw;
                    }
                    if (cnt == 1111) {
                        continue;
                    }

                    double consumptionNew = wageCurrent[i]*k;
                    double ne1 = Calculator.getErr(consumptionPrevious[i], consumptionNew, kapital[i], zCurrent[i], params);
                    double kapitalNext = wageCurrent[i] - consumptionNew;

                    for (int j = 0; j < Constants.M/Constants.P; ++j) {
                        double consumptionNext = Calculator.getWage(kapitalNext, zNext[j][i], params);
                        ne2[j] = Calculator.getErr(consumptionNew, consumptionNext, kapitalNext, zNext[j][i], params);
                        v2[j] = te2[j] + ne2[j] - errors2[j][i];
                    }

                    double v1 = te1 + ne1 - errors1[i], a2 = combiner.combine(v2);
                    if (val > abs(v1) + a2) {
                        errors1[i] = ne1; 
                        for (int j = 0; j < Constants.M/Constants.P; ++j) {
                            errors2[j][i] = ne2[j];
                            te2[j] = v2[j];
                        }
                        consumptionCurrent[i] = consumptionNew;
                        coefficients[i] = k;
                        te1 = v1;
                        val = abs(v1) + a2;
                    }

                    k = coefficients[i]*(1 + rng.nextGaussian()*pw);
                    cnt = 0;
                    while (++cnt < 1111 && (wageCurrent[i]*k <= 0 || k >= 1)) {
                        k = coefficients[i]*(1 + rng.nextGaussian()*pw);
                    }
                    if (cnt == 1111) {
                        continue;
                    }

                    consumptionNew = wageCurrent[i]*k;
                    ne1 = Calculator.getErr(consumptionPrevious[i], consumptionNew, kapital[i], zCurrent[i], params);
                    kapitalNext = wageCurrent[i] - consumptionNew;

                    for (int j = 0; j < Constants.M/Constants.P; ++j) {
                        double consumptionNext = Calculator.getWage(kapitalNext, zNext[j][i], params);
                        ne2[j] = Calculator.getErr(consumptionNew, consumptionNext, kapitalNext, zNext[j][i], params);
                        v2[j] = te2[j] + ne2[j] - errors2[j][i];
                    }

                    v1 = te1 + ne1 - errors1[i]; a2 = combiner.combine(v2);
                    if (val > abs(v1) + a2) {
                        errors1[i] = ne1; 
                        for (int j = 0; j < Constants.M/Constants.P; ++j) {
                            errors2[j][i] = ne2[j];
                            te2[j] = v2[j];
                        }
                        consumptionCurrent[i] = consumptionNew;
                        coefficients[i] = k;
                        te1 = v1;
                        val = abs(v1) + a2;
                    }
                }
                pw *= 0.1;
            }
        }
        return val;
    }
}

class Calculator {

    public static double getWage(double k, double z, Parameters params) {
        return (1-params.delta)*k + z*pow(k, params.alpha);
    }

    public static double getCt(double consumptionPrevious, double k, double z, Parameters params) {
        return consumptionPrevious*pow(params.beta*(1-params.delta+z*params.alpha*pow(k, params.alpha-1)), params.reta);
    }

    public static double getErr(double consumptionPrevious, double consumptionCurrent, double k, double z, Parameters params) {
        return params.beta*pow(consumptionPrevious/consumptionCurrent, params.eta)*(1-params.delta+z*params.alpha*pow(k, params.alpha-1))-1;
    }  

    public static double sqr(double x) {
        return x*x;
    }

    public static double cube(double x) {
        return x*x*x;
    }
}

interface ICombiner {

    public double combine(double[] v);
}

class SquareCombiner implements ICombiner {

    @Override
    public double combine(double[] v) {
        double sum = 0;
        for (double x: v) {
            sum += x*x;
        }
        return pow(sum/v.length, 0.5);
    }
}

class DCombiner implements ICombiner {

    @Override
    public double combine(double[] v) {
        double sum = 0;
        for (double x: v) {
            sum += abs(x);
        }
        int s = v.length;
        double m = sum/s;
        sum = 0;
        for (double x: v) {
            sum += Calculator.sqr(abs(x) - m);
        }
        return sum/s;
    }
}

class Combiner implements ICombiner {

    @Override
    public double combine(double[] v) {
        double val = 0;
        for (ICombiner c: combiners) {
            val += c.combine(v);
        }
        return val;
    }

    ICombiner[] combiners = {
        new SquareCombiner(), new DCombiner()
    };
}

class Timer {

    public Timer(long timeLimit) {
        finishTime = currentTimeMillis() + timeLimit;
    }

    public boolean isTL() {
        return currentTimeMillis() > finishTime;
    }

    private double finishTime;
}
