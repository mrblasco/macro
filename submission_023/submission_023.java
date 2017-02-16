/**
 * V2
 * 
 * https://community.topcoder.com/longcontest/?module=ViewStandings&rd=16765
 * 
 * @author smarsoll
 *
 */
public class NationSave {

	// private final static MyLogger LOGGER = new MyLogger();
	// protected final static org.slf4j.Logger LOGGER = org.slf4j.LoggerFactory.getLogger(NationSave.class);
	/** 0.000001? */
	public static final double VERY_SMALL = 0.00002;
	public static final int NB_PARAMS = 6;

	protected int N;
	protected int iT;
	protected int T;

	protected double beta, alpha, delta;

	protected double[] fact;
	protected double[] resultP1;

	/**
	 * Î² (beta) - Future utility discount factor<br>
	 * Î· (eta) - Consumption utility parameter<br>
	 * Î± (alpha) - Productivity of capital<br>
	 * Î´ (delta) - Depreciation rate of capital<br>
	 * Ï (rho) - Degree of correlation with previous shocks<br>
	 * Ï (sigma) - Standard deviation of Îµt<br>
	 * N - Number of concurrent simulations<br>
	 * T - Total number of time periods<br>
	 * 
	 * @param beta
	 * @param eta
	 * @param alpha
	 * @param delta
	 * @param rho
	 * @param sigma
	 * @param N
	 * @param T
	 * @return
	 */
	public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
		// LOGGER.info("SetEconomyParameters");

		fact = new double[NB_PARAMS];
		resultP1 = new double[N];

		fact[0] = beta;
		fact[1] = eta;
		fact[2] = alpha;
		fact[3] = delta;
		fact[4] = rho;
		fact[5] = sigma;

		this.beta = beta;
		this.alpha = alpha;
		this.delta = delta;
		this.N = N;
		this.T = T;
		this.iT = 0;
		return 0;
	}

	/**
	 * Kt - Invested capital at current time for each of the N concurrent simulations<br>
	 * Zt - Shock at current time for each of the N concurrent simulations<br>
	 * 
	 * ConsumptionDecisionRule must return the consumption for this time period for each of the N concurrent
	 * simulations. Consumption at each time period for each simulation must be real, non-negative, and not greater than
	 * the current wage.
	 * 
	 * @param Kt
	 * @param Zt
	 * @return
	 */
	public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt) {
		// LOGGER.info("Kt[0]=" + Kt[0] + " , Zt[0]=" + Zt[0]);
		double[] result = new double[N];
		double[] W_t = new double[N];
		// double[] result = Kt;
		// fact[9] = iT / (double) T;
		for (int i = 0; i < N; i++) {

			// Compute wage
			W_t[i] = (1 - delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
			// LOGGER.info("Kt[i]=" + Kt[i] + " , W_t=" + W_t);

			// fact[0] = Kt[i];
			// fact[1] = Zt[i];
			// fact[2] = W_t[i];

			// result[i] = (W_t[i] * 0.34);

			// S[sc=640619,15461196 at 2000, 3x[C 0,35 [], C 0,256 [6], C -0,032 [4]]]
			// result[i] = W_t[i] * (0.35 + 0.256 * factValues[6] - 0.032 * factValues[4]);
			// S[sc=673175,23434734 at 2000, 4x[C 0,448 [6], C 0,35 [], C -0,192 [8, 8, 8], C -0,032 [4]]]
			// result[i] = W_t[i] * (0.35 + 0.448 * fact[6] - 0.192 * fact[8] * fact[8] * fact[8] - 0.032 * fact[4]);

			// S[sc=684511.51178919 at 5000, 6x[C 0.576 [6], C 0.35 [], C -0.256 [4, 5, 8], C -0.128 [8, 8, 8], C 0.032
			// [3], C -0.032 [4]]]
			// result[i] = W_t[i]
			// * (0.35 + 0.576 * fact[6] - 0.256 * fact[4] * fact[5] * fact[8] - 0.128 * fact[8] * fact[8] * fact[8] +
			// 0.032
			// * fact[3] - 0.032 * fact[4]);

			// S[sc=713249.61733519 at 17000, 10x[C 0.57042768 [6], C -0.42834052 [4, 5, 8], C 0.35726662 [], C
			// -0.19234839 [8, 8, 8], C 0.0618062 [3], C 0.03906074 [3, 3], C -0.0354908 [6, 8, 8], C -0.032 [4], C
			// 0.02872478 [3, 6, 6], C -0.01470091 [4, 8]]]
			// result[i] = W_t[i]
			// * (0.35726662 + 0.57042768 * fact[3] - 0.42834052 * fact[1] * fact[2] * fact[5] - 0.19234839 * fact[5] *
			// fact[5]
			// * fact[5] + 0.0618062 * fact[2] + 0.03906074 * fact[2] * fact[2] - 0.0354908 * fact[3] * fact[5] *
			// fact[5]
			// - 0.032 * fact[1] + 0.02872478 * fact[0] * fact[3] * fact[3] - 0.01470091 * fact[1] * fact[5]);
			result[i] = W_t[i]
					* ( //
					0.40334336 + 0.22997859 * fact[3] + 0.21837595 * fact[3] * fact[3] + 0.15973829 * fact[3] * fact[4] + -0.10161004
							* fact[1] * fact[5] + -0.020384 * fact[1] * fact[1] * fact[5]
			);
			// Keep only valid solutions
			result[i] = Math.max(result[i], VERY_SMALL * W_t[i]);
			// to have a valid result
			result[i] = Math.min(result[i], W_t[i] * (1.0 - VERY_SMALL));
		}

		// LOGGER.info("Kt[0]=" + Kt[0] + " , Zt[0]=" + Zt[0] + " , " + result[0] + " for " + iT);
		// int s = 0;
		// LOGGER.info("Kt[" + s + "]=" + Kt[s] + " , Zt[]=" + Zt[s] + " , W_t[]=" + W_t[s] + " , " + result[s] +
		// " for " + iT);
		iT++;
		return result;
	}


}

class MyLogger {
	public final void info(String string) {
		System.out.println(string);
	}

	public void warn(String string) {
		System.out.println(string);
	}

	public final boolean isInfoEnabled() {
		return false;
	}

	public final void error(String string) {
	}
}