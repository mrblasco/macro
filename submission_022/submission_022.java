public class NationSave {
    //  economy parameters
        	double beta = 0;
            double eta = 1;
            double alpha = 2;
            double delta = 3;
            double rho = 4;
            double sigma = 5;
  int stopprint=1;
    
    public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
    {
        double[] ans= new double[Kt.length];
        //System.out.println(Kt.length); 
        
        double W_t,factor;
        
        for(int i=0;i<Kt.length;i++)
        {
            W_t = (1 - delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
            
          factor = 1.662354- 0.765820 * beta - 0.118725 *eta - 1.074698 * alpha + 0.483890* delta - 0.335990 *sigma;
           
           if(factor<=0) factor =0.1;
           if(factor>=1) factor= 0.9;
            ans[i] = W_t*factor;
        }
        return ans;
    }
    
    public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T)
    {
    this.alpha =alpha;
    this.beta = beta;
    this.delta = delta;
    this.eta = eta;
    this.rho = rho;
    this.sigma =sigma;
    return 0;
    }
}