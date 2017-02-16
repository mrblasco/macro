import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

public class NationSave {
    
    int N = 1;
    double beta = 0.0;
    double eta = 0.0;
    double alpha = 0.0;
    double delta = 0.0;
    double rho = 0.0;
    double sigma = 0.0;
    int T = 1;
    double[] Wt;
    
    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt) {
        double [] ret = new double[N];
        computeWt(Kt, Zt);
        
        for(int i = 0; i<N;i++) {
            double cCandidate = 0.0;
            double cMin = 0.0;
            double cMax = Wt[i];
            if (Wt[i] <= 0.0) {
                ret[i] = 0.0;
                continue;
            }
            ret [i] = Wt[i]/2.0;
        }
        return ret;
    }

    private void computeWt(double[] Kt, double[] Zt) {
        Wt = new double[N];        
        for(int i = 0; i<N;i++) {
            Wt[i] = (1-delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
        }
    }
    
    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
        this.N = N;
        this.beta = beta;
        this.eta = eta;
        this.alpha = alpha;
        this.delta = delta;
        this.rho = rho;
        this.sigma = sigma;
        this.T = T;
        return 0;
    }
    
    public static void main(String[] args) throws NumberFormatException, IOException {
        NationSave nationSave = new NationSave();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        PrintWriter writer = new PrintWriter(System.out);
        
        int X = Integer.parseInt(reader.readLine());
        for (int xx = 0; xx < X; xx++ ){
            double beta = Double.parseDouble(reader.readLine());
            double eta = Double.parseDouble(reader.readLine());
            double alpha = Double.parseDouble(reader.readLine());
            double delta = Double.parseDouble(reader.readLine());
            double rho = Double.parseDouble(reader.readLine());
            double sigma = Double.parseDouble(reader.readLine());
            int N = Integer.parseInt(reader.readLine());
            int T = Integer.parseInt(reader.readLine());
            
            double [] Kt = new double[N];
            double [] Zt = new double[N];
            
            nationSave.SetEconomyParameters(beta, eta, alpha, delta, rho, sigma, N, T);
            
            for (int tt = 0; tt<T; tt++) {
                for (int i = 0; i<N; i++)  {
                    Kt[i] = Double.parseDouble(reader.readLine());
                }
                
                for (int i = 0; i<N; i++)  {
                    Zt[i] = Double.parseDouble(reader.readLine());
                }
                
                double [] Ct = nationSave.ConsumptionDecisionRule(Kt, Zt);
                
                for (int i = 0; i<N; i++)  {
                    writer.println(Ct[i]);
                }
                writer.flush();            
            }
        }       
    }

}
