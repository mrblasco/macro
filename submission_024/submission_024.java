public class NationSave 
{
	public int stepcnt=0;
	public double beta,eta,alpha,delta,rho,sigma; 
	public int N,T;
	public double det_factor=0.001;
	public double multiplier=0.4;
	// Consumption
	public double[][] C;
				// Capital
	public double[][] K;
				// Shock
	public double[][] Z;
	
	public double[][] W;
	public double[][] multiplier_matrix;
				
	public double scoreAnswer(int Tt,double[] Ct,double[] Zt,double[] Kt)
	{
		// Scoring
		double eer = 0;
		for(int t = 0; t < 8; t++)
		{
			eer+= Math.abs( beta * Math.pow(Ct[t]/ Ct[t + 1], eta) * (1 - delta + Zt[t + 1] * alpha * Math.pow(Kt[t + 1], alpha - 1)) - 1);
			
		}
		return 1000000/(1 + (double)eer/8);
		
		
	}
	
	
	public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt)
	{
		
	
		double[] Ct=new double[N];
		double Wt=0;
		if(stepcnt%10==0 && stepcnt>0)
		{
			double tmp_multiplier=0.05;
			double tmpscore=0;
			int Ttmp=stepcnt-(stepcnt/10)*10;
			double[] Ctmp=new double[10];
			double[] Ztmp=new double[10];
			double[] Ktmp=new double[10];
			for(int i=0;i<N;i++)
			{
				tmp_multiplier=0.05;
				for(int j=0;j<15;j++)
				{					
					tmpscore=0;
					Ktmp[0]=K[(stepcnt/10-1)*10][i];
					Ztmp[0]=Z[(stepcnt/10-1)*10][i];
					int m=0;
					for(int k=(stepcnt/10-1)*10;k<stepcnt-1;k++)
					{
						Ztmp[m]=Z[k][i];
						double Wtmp=(1-delta)*Ktmp[m]+Z[k][i]*Math.pow(Ktmp[m], alpha);
						Ctmp[m]=Wtmp*tmp_multiplier;
						Ktmp[m + 1] = Wtmp - Ctmp[m];
						m++;
					}
					tmpscore=scoreAnswer(Ttmp, Ctmp, Ztmp, Ktmp);
					if(tmpscore>multiplier_matrix[i][1])
					{
						multiplier_matrix[i][1]=tmpscore;
						multiplier_matrix[i][0]=tmp_multiplier;
					}
					tmp_multiplier=tmp_multiplier+0.05;
				}
				
			}
			
		}
		for(int i=0;i<N;i++)
		{
			
			Wt=(1-delta)*Kt[i]+Zt[i]*Math.pow(Kt[i], alpha);
			
			Ct[i]=Wt*multiplier_matrix[i][0];
			C[stepcnt][i]=Ct[i];
			W[stepcnt][i]=Wt;
			K[stepcnt][i]=Kt[i];
			Z[stepcnt][i]=Zt[i];
				//multiplier=multiplier-det_factor;
			
		}
		stepcnt++;
		return Ct;
		
	}
	
	public int SetEconomyParameters(double beta1, double eta1, double alpha1, double delta1, double rho1, double sigma1, int N1, int T1)
	{
		N=N1;
		T=T1;
		beta=beta1;
		eta=eta1;
		alpha=alpha1;
		delta=delta1;
		rho=rho1;
		sigma=sigma1;
		// Consumption
		C = new double[T + 1][N];
					// Capital
		K = new double[T + 1][N];
					// Shock
		Z = new double[T + 1][N];
		
		W = new double[T + 1][N];
		multiplier_matrix=new double[N][2];
		det_factor=0.5/N;
		multiplier=0.7;
		if(sigma<0.25 || eta<1.2)
		{
			for(int i=0;i<N;i++)
				multiplier_matrix[i][0]=0.7;
		}
		else if(alpha<.18)
		{
			for(int i=0;i<N;i++)
				multiplier_matrix[i][0]=0.25;
		}
		else
		{
			for(int i=0;i<N;i++)
				multiplier_matrix[i][0]=0.1;
		}
		
		return 0;
	}
	
}
