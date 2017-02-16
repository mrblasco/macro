// package org.onja.topcoder.nationsave;

/**
 * NationSave
 * */
public class NationSave {

	public double beta; // Future utility discount factor
	public double eta; // Consumption utility parameter
	public double alpha; // Productivity of capital
	public double delta; // Depreciation rate of capital
	public double rho; // Degree of correlation with previous shocks
	public double sigma; // Standard deviation of ÃÂµt
	public int N; // Number of concurrent simulations
	public int T; // Total number of time periods
	
	public int t = 0; // Current Time. 
	public double W[][]; // Wages
	
	/**
	 * SetEconomyParameters
	 * @param double beta Future utility discount factor
	 * @param double eta Consumption utility parameter
	 * @param double alpha Productivity of capital
	 * @param double delta Depreciation rate of capital
	 * @param double rho Degree of correlation with previous shocks
	 * @param double sigma Standard deviation of ÃÂµt
	 * @param int N Number of concurrent simulations
	 * @param int T Total number of time periods
	 * @return int Ignored. 
	 * */
	public int SetEconomyParameters(double beta, double eta, double alpha, double delta, 
			double rho, double sigma, 
			int N, int T) {		
		this.beta = beta;
		this.eta = eta;
		this.alpha = alpha;
		this.delta = delta;
		this.rho = rho;
		this.sigma = sigma;
		this.N = N;
		this.T = T;
		
		this.W = new double[N][T];
		return 0;
	}

	/**
	 * ConsumptionDecisionRule
	 * @param double[] k Invested capital at current time for each of the N concurrent simulations
	 * @param double[] z Shock at current time for each of the N concurrent simulations
	 * */
	public double[] ConsumptionDecisionRule(double[] k, double[] z) {
		double[] result = new double[N];
		for (int i = 0; i < N; i++) {
			W[i][t] = (1 - delta) * k[i] + z[i] * Math.pow(k[i], alpha); // Compute wage
			result[i] = W[i][t] * 0.25;
		}

		t++;
		return result;
	}
	
}