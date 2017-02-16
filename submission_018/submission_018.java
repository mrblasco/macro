import java.util.*;

public class NationSave{
  double beta, eta, alpha, delta, rho, sigma;
  int N, T;
  Test test;
  
  public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt){
    double[] Ct = new double[N];
    
    for (int n = 0; n < N; n++){
      double Wt = (1 - delta)*Kt[n] + Zt[n]*Math.pow(Kt[n], alpha);
      Ct[n] = Wt*test.bestFrac;
    }
    
    return Ct;
  }
  
  public int SetEconomyParameters(double beta, double eta, double alpha, 
                                  double delta, double rho, double sigma, 
                                  int N, int T){
    this.beta = beta;   // future utility discount factor
    this.eta = eta;     // consumption utility parameter
    this.alpha = alpha; // productivity of capital
    this.delta = delta; // depreciation rate of capital
    this.rho = rho;     // degree of correlation with previous shocks
    this.sigma = sigma; // standard deviation of epsilon_t
    this.N = N;         // number of concurrent simulation
    this.T = T;         // total number of time periods
    test = new Test();
    test.run();
    return 0;
  }

  class Test{
    final static double minW = 0.0001;
    double bestScore, bestFrac;
    
    double[][] shocks = new double[T][N];
    
    double K0 = 1, Z0 = 1;
    double[][] C = new double[T + 1][N];
    double[][] K = new double[T + 1][N];
    double[][] Z = new double[T + 1][N];
    double[] Wt = new double[N];
    
    void run(){
      Random rng = new Random(1);
      for(int t = 0; t < T; t++){
        for(int n = 0; n < N; n++){
          shocks[t][n] = rng.nextGaussian()*sigma;
        }
      }
      
      bestFrac = 0.5;
      bestScore = testFrac(bestFrac);
      testRange(0, 1);
    }
    
    void testRange(double x, double w){
      if (w < minW) return;
      double f = x + w/4, score = testFrac(f);
      if (score > bestScore){
        bestScore = score; bestFrac = f; testRange(x, w/2);
      }else{
        f+= w/2; score = testFrac(f);
        if (score > bestScore){
          bestScore = score; bestFrac = f; testRange(x+w/2, w/2);
        }else
          testRange(x+w/4, w/2);
      }
    }
    
    public double testFrac(double frac){
      Arrays.fill(K[0], K0);
      Arrays.fill(Z[0], Z0);
			
      for(int t = 0; t <= T; t++){
        for(int n = 0; n < N; n++){
          if (t != 0)
            Z[t][n] = Math.exp(rho * Math.log(Z[t - 1][n]) + shocks[t - 1][n]);
          Wt[n] = (1 - delta) * K[t][n] + Z[t][n] * Math.pow(K[t][n], alpha);
        }
        if (t != T){
          for (int n = 0; n < N; n++){ 
            C[t][n] = frac*Wt[n];
            K[t + 1][n] = Wt[n] - C[t][n];
          }
        }else{
          C[t] = Wt;
        }
      }
      
      // Scoring
      double eer = 0;
      for(int t = 0; t < T; t++){
        double eer_t = 0;
        for(int i = 0; i < N; i++){
          double eer_t_i = beta*Math.pow(C[t][i]/C[t+1][i], eta)*
               (1 - delta + Z[t+1][i]*alpha*Math.pow(K[t+1][i], alpha-1)) - 1;
          eer_t+= eer_t_i;
        }
        eer+= Math.abs(eer_t/N);
      }
      return 1000000.0/(1 + eer/T);
    }
  }
  
}
