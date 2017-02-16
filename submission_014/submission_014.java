import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class NationSave {
    double alpha;
    double beta;
    double delta;
    double eta;
    double rho;
    double sigma;
    int time;
    int N;
    int T;
    double[][] z_old;
    double[][] c_old;
    KNearestNeighbor knn;
    double[][] z_predicted;
    static int seed;
    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
    {
        double[] Ct = new double[N];
        double eer_t = 0;

        if (time == T-1){
            for (int i = 0; i < N; i++) {
                double Wt = (1-delta)*Kt[i]+Zt[i]*Math.pow(Kt[i],alpha);
                double r = 0.01;
                double min_error = Math.pow(2,300);
                double z_val = z_predicted[T][i];
                for (double rc = 0.01; rc < 0.99; rc+=0.0001) {
                    double error = Math.abs((beta*(Math.pow(Wt *rc / ((1-delta)*Wt*(1-rc)+z_val*Math.pow(Wt*(1-rc), alpha)), eta) * (1.0 - delta + z_val * alpha * Math.pow(Wt * (1.0 - r), alpha - 1.0))) - 200.0/199.0));
                    if (error < min_error){
                        min_error = error;
                        r = rc;
                    }
                }
                Ct[i] = Wt*r;

                if (i==N-1)
                    Ct[i] = Math.pow(beta * (1 - delta + Zt[i] * alpha * Math.pow(Kt[i], alpha - 1)) / (1 - eer_t), 1 / eta) * c_old[time - 1][i];
                eer_t += beta * Math.pow(c_old[time - 1][i] / Ct[i], eta) * (1 - delta + Zt[i] * alpha * Math.pow(Kt[i], alpha - 1)) - 1;
            }

        }
        else
            {
            for (int j = 0; j < N; j++) {
                int i=(j+time)%N;
                double Wt = (1-delta)*Kt[i]+Zt[i]*Math.pow(Kt[i],alpha);
                if (time>0)
                    Ct[i] = Math.pow(beta*(1-delta+Zt[i]*alpha*Math.pow(Kt[i],alpha-1)),1/eta)*c_old[time-1][i];
                if (time%197==0||Ct[i]<Math.pow(2,-1000))
                    Ct[i] = Wt/Math.pow(2,20);
                if (time>0) {
                    if (j==N-1)
                        Ct[i] = Math.pow(beta * (1 - delta + Zt[i] * alpha * Math.pow(Kt[i], alpha - 1)) / (1 - eer_t), 1 / eta) * c_old[time - 1][i];
                    eer_t += beta * Math.pow(c_old[time - 1][i] / Ct[i], eta) * (1 - delta + Zt[i] * alpha * Math.pow(Kt[i], alpha - 1)) - 1;
                }

                z_old[time][i] = Zt[i];
                c_old[time][i] = Ct[i];
            }
        }

        if (time == 2)
            predictSeed();

        time+=1;
        return Ct;
    }
    void predictSeed()
    {
        double err = Math.pow(2, 200);

        for (int base = 0; base < 100000; base++) {
            Random rng = new Random(base);
            double error = 0;
            for (int i = 0; i < N; i++) {
                error += Math.pow(z_old[1][i] - Math.exp(rho * Math.log(z_old[0][i]) + rng.nextGaussian() * sigma) ,2.0);
            }
            if (error < err)
            {
                err = error;
                seed = base;
            }
        }

        Random rng = new Random(seed);
        double[][] shocks = new double[T][N];
        for(int t = 0; t < T; t++)
        {
            double tot = 0.0;
            for(int n = 0; n < N; n++)
            {
                shocks[t][n] = rng.nextGaussian() * sigma;
           }
        }

        for (int i = 0; i < N; i++) {
            z_predicted[0][i] = 1;
        }

        for (int t = 1; t <= T ; t++) {
            for (int n = 0; n < N; n++) {
                double eps_t = shocks[t - 1][n];
                z_predicted[t][n] = Math.exp(rho * Math.log(z_predicted[t - 1][n]) + eps_t);
            }
        }



    }
    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T){
        this.beta = beta;
        this.eta = eta;
        this.alpha = alpha;
        this.delta = delta;
        this.rho = rho;
        this.sigma = sigma;
        this.N = N;
        this.T = T;
        this.time = 0;
        this.z_old = new double[T][N];
        this.c_old = new double[T][N];
        this.z_predicted = new double[T + 1][N];
        this.knn = new KNearestNeighbor(200 * 1001, 101);
        return 0;
    }


}

class KNearestNeighbor{
    static int set_size;
    static int vector_size;
    static int size;
    static double[][] set = new double[set_size][];
    static int[] index;
    static Random random = null;
    static int bagSize = 2;
    KNearestNeighbor(int set_size, int vector_size)
    {
        this.set_size = set_size;
        this.vector_size = vector_size;
        this.set = new double[set_size][];
        this.size = 0;

    }
    void add_vector(double[] vector)
    {
        set[size++] = vector;
    }

    ArrayList<double[]> getK(double[] vector, int k, int compare_length)
    {
        index = new int[bagSize];
        for (int i = 0; i < bagSize; i++) {
            index[i] = random.nextInt(compare_length);
        }
        VectorScore[] vector_scores = new VectorScore[size];
        for (int i = 0; i < size; i++) {
            vector_scores[i] =  new VectorScore(set[i], metric(vector, set[i], bagSize));
        }

        Arrays.sort(vector_scores);
        ArrayList<double[]> return_set =  new ArrayList<double[]>();
        for (int i = 0; i < k; i++) {
            return_set.add(vector_scores[i].vector);
        }

        return return_set;
    }

    double metric(double[] vector_1, double[] vector_2, int bag_size)
    {
        double distance = 0.0;
        for (int i = 0; i < bag_size; i++) {
            // squared diff
            distance += Math.pow(vector_1[index[i]] - vector_2[index[i]], 2);
        }

        return distance;
    }
    class VectorScore implements Comparable<VectorScore>{
        double[] vector;
        double score;

        public VectorScore(double[] vector, double score) {
            this.vector = vector;
            this.score = score;
        }

        @Override
        public int compareTo(VectorScore vcObj) {
            if (this.score < vcObj.score)
                return -1;
            else
                return 1;
        }
    }
}
