import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.*;
import java.util.Queue;
import java.util.Random;
import java.util.regex.*;
import java.io.*;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class NationSave {
    double beta; double eta; double alpha; double delta; double rho; double sigma; int N; int T;
    
    double steadystate = (1/beta - 1 + delta*(1-alpha))/alpha;
    double ssfrac = steadystate/(1+steadystate);
    double spendpercent = 50.0/100;
    double last[] = new double[N];
    int count = 0;
    double[] ConsumptionDecisionRule(double[] Kt, double[] Zt){
        count++;
        //if (count<=2)System.err.println(Zt[1] + " " + Kt[1]);
        //if (count<10)System.err.println(count + "   -   " + Zt[0]);
        double lower = Math.exp(rho*Math.log(Zt[0]) - 1.645*sigma);
        double upper = Math.exp(rho*Math.log(Zt[0]) + 1.645*sigma);
        //if (count<10)System.err.println(lower + "   -   " + upper);
        double W_t[] = new double[Kt.length];
        double Ct[] = new double[Kt.length];
        for (int n=0;n<W_t.length;n++){
	W_t[n] = (1 - delta) * Kt[n] + Zt[n] * Math.pow(Kt[n], alpha);
        Ct[n] = spendpercent*W_t[n];
       
        double steadystate = (1/beta - 1 + delta*(1-alpha))/alpha;
        double ssfrac = steadystate/(1+steadystate);

        W_t[n] = (1 - delta) * Kt[n] + Zt[n] * Math.pow(Kt[n], alpha);
        //Ct[n] = ssfrac*W_t[n];
        
        //xt is log Zt
        //tau is eta
        double tau = eta;
        double phi = 1 + 1/beta + ((1-alpha)/tau)*(1-(1-delta)*beta)*steadystate;
        
        //solutions to x^2 - phi x + 1/beta
        //x = -b +- sqrt(b^2 - 4ac) / 2a    a = 1, b = -phi  c = 1/beta
        
        double xplus = (phi + Math.sqrt(phi*phi - 4/beta))/2.0;
        double xminus = (phi - Math.sqrt(phi*phi - 4/beta))/2.0;
        
        double lambda = 0;
        if (Math.abs(xplus)<1)lambda = xplus;
        else lambda = xminus;
        //System.err.println(" - " + lambda);
        double Kstar = Math.pow(beta*alpha/(1-(1-delta)*beta), 1/(1-alpha));
        double q = beta*((1-rho)*(steadystate + delta) + rho*beta*(1/beta + delta - 1) * steadystate / tau) * Kstar;
        //if (count>1){
        
            double xt = Math.log(Zt[n]);
            double K_t = (1-lambda)*Kstar + lambda*Kt[n] + q*lambda * Math.log(Zt[n]) /(1-beta*rho*lambda);
            
            //double K_t = Math.pow(Kstar, 1-lambda)*Math.pow(Kt[n], lambda);
            
            //W_t[n] = (1 - delta) * K[t][n] + Z[t][n] * Math.pow(K[t][n], alpha);
            
            //System.err.println(K_t);
            //if (K_t < 0 ) K_t = 0;
            //if (K_t > W_t[n] ) K_t = W_t[n];
            double C_t = W_t[n] - K_t;
            Ct[n] = C_t;
            if (K_t < 0 || K_t > W_t[n])Ct[n] = ssfrac * W_t[n];
            
            if (Kt[0]!=1){
            double opt = last[n]*Math.pow((beta*(1-delta+Zt[n]*alpha*Math.pow(Kt[n],alpha-1))),1/eta);
            Ct[n] = opt;
            if (Ct[n] > W_t[n] || opt/last[n] < 0.01 || opt < 1e-320){
                //if(sigma < 0.6)
                    Ct[n] = ssfrac*W_t[n];
                    Ct[n] = init;
                    if (Ct[n] > W_t[n])Ct[n] = bestval * W_t[n];
          	//else Ct[n] = 0.5 * W_t[n];
                //Ct[n] = 0.01*W_t[n];
            }
            //if (test==4)System.err.println(opt);
            }
            else {
                if(sigma < 0.6)Ct[n] = ssfrac*W_t[n];
          	else Ct[n] = 0.5 * W_t[n];
                Ct[n] = bestval * W_t[n];
                init = Ct[n];
            }
        }
        
        
        
        last = Ct;
        return Ct;
    }
    double bestval = 0.01;
    double init = 0;
    double scoreAnswer(int N, int T, double beta, double eta, double alpha, double delta, double[][] C, double[][] K, double[][] Z)
    {
            // Scoring
            double eer = 0;
            for(int t = 0; t < T; t++)
            {
                    double eer_t = 0;
                    for(int i = 0; i < N; i++)
                    {
                            double eer_t_i = beta * Math.pow(C[t][i] / C[t + 1][i], eta) * (1 - delta + Z[t + 1][i] * alpha * Math.pow(K[t + 1][i], alpha - 1)) - 1;
                            eer_t += eer_t_i;
                    }
                    eer_t /= N;
                    eer += Math.abs(eer_t);
            }
            eer = 1 / (1 + eer / T);

            double score = eer * 1000000;

            return score;
    }
    double simulate(){
        final int N = 200;
        final int T = 2000;
        final double K0 = 1;
        final double Z0 = 1; //get these
        
        int numsimu = 10;
        double steadystate = (1/beta - 1 + delta*(1-alpha))/alpha;
        double ssfrac = steadystate/(1+steadystate);
        double[] vals = {1e-100, 1e-50, 1e-37, 1e-25, 1e-16, 0.000000001, 0.000001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.7};
        double bestscore = 0;
        double best = 0.01;
        double[] scores = new double[vals.length];
        for (int sim = 0; sim<numsimu; sim++){
        for (int iter = 0;iter<vals.length;iter++){
            count2 = 0;
        	// Shock random numbers
        	double[][] shocks = new double[T][N];
        	Random rng = new Random(10*sim + iter + 1);
			for(int t = 0; t < T; t++)
			{
				for(int n = 0; n < N; n++)
				{
					shocks[t][n] = rng.nextGaussian() * sigma;
				}
			}
			
			// Consumption
			double[][] C = new double[T + 1][N];
			// Capital
			double[][] K = new double[T + 1][N];
			// Shock
			double[][] Z = new double[T + 1][N];
			
			// Set initial capital
			Arrays.fill(K[0], K0);
			// Set initial shock
			Arrays.fill(Z[0], Z0);
			
			for(int t = 0; t <= T; t++)
			{
				// Wage
				double[] W_t = new double[N];
				for(int n = 0; n < N; n++)
				{
					// Update shock
					if(t != 0) //Already set for t = 0
					{
						double eps_t = shocks[t - 1][n];
						Z[t][n] = Math.exp(rho * Math.log(Z[t - 1][n]) + eps_t);
					}
					// Compute wage
					W_t[n] = (1 - delta) * K[t][n] + Z[t][n] * Math.pow(K[t][n], alpha);
				}
				if(t != T)
				{
					/*for(int i = 0; i < N; i++)
					{
						writer.println(K[t][i]);
					}
					for(int i = 0; i < N; i++)
					{
						writer.println(Z[t][i]);
					}
					
					writer.flush();
					
					for(int i = 0; i < N; i++)
					{
						C[t][i] = Double.parseDouble(reader.readLine());
					}
					*/
                                        C[t] = this.ConsumptionDecisionRuleParam(K[t], Z[t], vals[iter]);
					for(int i = 0; i < N; i++)
					{
                                            /*
						//Sanity check C[t][i]
						if(C[t][i] < 0)
						{
							printMessage("Consumption cannont be negative.");
							printMessage("Score = 0");
							return;
						}
						if(C[t][i] > W_t[i])
						{
							printMessage("Consumption cannont exceed wage.");
							printMessage("Score = 0");
							return;
						}
						if(Double.isNaN(C[t][i]))
						{
							printMessage("Consumption must be a real number.");
							printMessage("Score = 0");
							return;
						}
						*/
						//Update capital for next period
						K[t + 1][i] = W_t[i] - C[t][i];
					}
				}
				else
				{
					//On last time period, consume all of the wage
					C[t] = W_t;
				}
			}
			
			// Scoring
        double testScore = scoreAnswer(N, T, beta, eta, alpha, delta, C, K, Z);
        scores[iter] += testScore;
        if (bestscore<scores[iter]){
            bestscore = scores[iter];
            best = vals[iter];
        }
        }
        }
        for (int iter = 0;iter < vals.length; iter ++){
            System.err.println(iter + " - " + scores[iter]/numsimu);
        }
        return best;
    }
    
    double initval = 0;
    int count2 = 0;
    double[] ConsumptionDecisionRuleParam(double[] Kt, double[] Zt, double initval){
        
        count++;
        //if (count<=2)System.err.println(Zt[1] + " " + Kt[1]);
        //if (count<10)System.err.println(count + "   -   " + Zt[0]);
        double lower = Math.exp(rho*Math.log(Zt[0]) - 1.645*sigma);
        double upper = Math.exp(rho*Math.log(Zt[0]) + 1.645*sigma);
        //if (count<10)System.err.println(lower + "   -   " + upper);
        double W_t[] = new double[Kt.length];
        double Ct[] = new double[Kt.length];
        for (int n=0;n<W_t.length;n++){
	W_t[n] = (1 - delta) * Kt[n] + Zt[n] * Math.pow(Kt[n], alpha);
        Ct[n] = spendpercent*W_t[n];
       
        double steadystate = (1/beta - 1 + delta*(1-alpha))/alpha;
        double ssfrac = steadystate/(1+steadystate);

        W_t[n] = (1 - delta) * Kt[n] + Zt[n] * Math.pow(Kt[n], alpha);
        //Ct[n] = ssfrac*W_t[n];
        
        //xt is log Zt
        //tau is eta
        double tau = eta;
        double phi = 1 + 1/beta + ((1-alpha)/tau)*(1-(1-delta)*beta)*steadystate;
        
        //solutions to x^2 - phi x + 1/beta
        //x = -b +- sqrt(b^2 - 4ac) / 2a    a = 1, b = -phi  c = 1/beta
        
        double xplus = (phi + Math.sqrt(phi*phi - 4/beta))/2.0;
        double xminus = (phi - Math.sqrt(phi*phi - 4/beta))/2.0;
        
        double lambda = 0;
        if (Math.abs(xplus)<1)lambda = xplus;
        else lambda = xminus;
        //System.err.println(" - " + lambda);
        double Kstar = Math.pow(beta*alpha/(1-(1-delta)*beta), 1/(1-alpha));
        double q = beta*((1-rho)*(steadystate + delta) + rho*beta*(1/beta + delta - 1) * steadystate / tau) * Kstar;
        //if (count>1){
        
            double xt = Math.log(Zt[n]);
            double K_t = (1-lambda)*Kstar + lambda*Kt[n] + q*lambda * Math.log(Zt[n]) /(1-beta*rho*lambda);
            
            //double K_t = Math.pow(Kstar, 1-lambda)*Math.pow(Kt[n], lambda);
            
            //W_t[n] = (1 - delta) * K[t][n] + Z[t][n] * Math.pow(K[t][n], alpha);
            
            //System.err.println(K_t);
            //if (K_t < 0 ) K_t = 0;
            //if (K_t > W_t[n] ) K_t = W_t[n];
            double C_t = W_t[n] - K_t;
            Ct[n] = C_t;
            if (K_t < 0 || K_t > W_t[n])Ct[n] = ssfrac * W_t[n];
            
            if (Kt[0]!=1){
            double opt = last[n]*Math.pow((beta*(1-delta+Zt[n]*alpha*Math.pow(Kt[n],alpha-1))),1/eta);
            Ct[n] = opt;
            if (Ct[n] > W_t[n] || opt/last[n] < 0.01 || opt < 1e-320){
                //if(sigma < 0.6)
                    Ct[n] = ssfrac*W_t[n];
                    Ct[n] = init;
                    if (Ct[n] > W_t[n])Ct[n] = bestval * W_t[n];
          	//else Ct[n] = 0.5 * W_t[n];
                //Ct[n] = 0.01*W_t[n];
            }
            //if (test==4)System.err.println(opt);
            }
            else {
                if(sigma < 0.6)Ct[n] = ssfrac*W_t[n];
          	else Ct[n] = 0.5 * W_t[n];
                Ct[n] = initval * W_t[n];
                init = Ct[n];
            }
        }
        
        
        
        last = Ct;
        return Ct;
    }
    
    int test = 0;
    int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T){
        count = 0;
        test++;
        System.err.println("Testcase - " + test);
        //System.err.println(beta + " , " + eta + " , " + alpha);
        this.beta = beta;
        this.eta = eta;
        this.alpha = alpha;
        this.delta = delta;
        this.rho = rho;
        this.sigma = sigma;
        this.N = N;
        this.T = T;
        

        bestval = simulate();
        System.err.println("value chosen = " + bestval);
        return 0;
    }
            
            
    public static void main(String argsv[]){
        NationSave ns = new NationSave();
        Scanner scan = new Scanner(System.in);
         int X = scan.nextInt();
        for (int w=0;w<X;w++){
           double beta = scan.nextDouble();
           double eta = scan.nextDouble();
           double alpha = scan.nextDouble();
           double delta = scan.nextDouble();
           double rho = scan.nextDouble();
           double sigma = scan.nextDouble();
           int N = scan.nextInt();
           int T = scan.nextInt();
           ns.SetEconomyParameters(beta, eta, alpha, delta, rho, sigma, N, T);

           for (int t=0;t<T;t++){
               double[] Kt = new double[N];
               double[] Zt = new double[N];

              for (int i=0;i<N;i++)
                 Kt[i] = scan.nextDouble();
              for (int i=0;i<N;i++)
                 Zt[i] = scan.nextDouble();

              double Ct[] = ns.ConsumptionDecisionRule(Kt, Zt);

              for (int i=0;i<N;i++)
                 System.out.println(Ct[i]);
            }
        }
    }
}