import java.util.Random;

// How Much Can a Nation Save - coolsmurf
public class NationSave {

	private static final int NUMBER_OF_SA = 3;
	private static final long TIME_PER_SA = 30 * 1000;
	private static final int ROLLS = 10000;

	private int N, T, t;
	private double beta, eta, alpha, delta, rho, sigma;
	private double[] last;

	private static final int PERIOD = 5;
	private static final double SHRINK = 300;

	public double[] ConsumptionDecisionRule(double[] Kt, double[] Zt) {

		double[] result = new double[N];
		if (t == 0) {
			for (int i = 0; i < N; i++) {
				double w = (1 - delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
				result[i] = w * Math.pow(10, -SHRINK);
			}
		} else if (t == T - 1) {
			for (int i = 0; i < N; i++) {
				result[i] = chooseLastLegal(Kt[i], Zt[i]);
			}
			result = simulatedAnnealing(Kt, Zt, result);
			result[0] = last[0] * Math.pow(beta * (1 - delta + alpha * Zt[0] * Math.pow(Kt[0], alpha - 1)) / N, 1 / eta);

		} else {
			result = smooth(Kt, Zt);
		}
		// Catastrophy recovery
		for (int i = 0; i < N; i++) {
			double w = (1 - delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
			if (result[i] > 0.99 * w) {
				result[i] = 0.99 * w;
			}
		}
		t++;
		last = result;
		return result;
	}

	private double[] smooth(double[] Kt, double[] Zt) {

		double meanScale = 0.0;
		for (int n = 0; n < N; n++) {
			meanScale -= Math.log10(last[n]) / N;
		}
		double multiplier = Math.pow(10, meanScale);
		if (t >= T - 2 || meanScale <= 5) {
			System.out.println("At t = " + t + ", mean scale is " + meanScale);
		}

		double[] result = new double[N];
		double[] wages = new double[N];
		double sum = 0.0;
		double sumLow = 0.0;
		double sumHigh = 0.0;
		for (int i = 0; i < N; i++) {
			wages[i] = (1 - delta) * Kt[i] + Zt[i] * Math.pow(Kt[i], alpha);
			double slope = beta * Math.pow(last[i] * multiplier / wages[i], eta) * (1 - delta + alpha * Zt[i] * Math.pow(Kt[i], alpha - 1)) / N;
			sum += slope;
			if ((t % PERIOD == i % PERIOD)) {
				sumLow += slope;
			} else {
				sumHigh += slope;
			}
		}
		if (eta * (SHRINK - meanScale) * Math.log(10) + Math.log(sum) > 0.0) {
			for (int i = 0; i < N; i++) {
				result[i] = wages[i] * Math.pow(10, -meanScale) * Math.pow(sum, 1.0 / eta);
			}
		} else {
			double B = Math.pow(sumHigh / (1.0 - sumLow * Math.pow(10, eta * (SHRINK - meanScale))), 1.0 / eta);
			for (int i = 0; i < N; i++) {
				if ((t % PERIOD == i % PERIOD)) {
					result[i] = wages[i] * Math.pow(10, -SHRINK);
				} else {
					result[i] = wages[i] * Math.pow(10, -meanScale) * B;
				}
			}
		}
		return result;
	}

	private double chooseLastLegal(double k, double z) {

		double w = (1 - delta) * k + z * Math.pow(k, alpha);
		int steps = 100;
		double bestPhiError = Double.MAX_VALUE;
		double bestS = 0.3;
		double goal = N / (N - 1.0);
		for (int i = 5; i < steps - 5; i++) {
			double s = 1.0 * i / steps;
			double p = 0.0;
			int S = 100;
			double eps = 3.0 / S;
			for (int sim = -S; sim <= S; sim++) {
				double x = eps * sim;
				double wp = eps * Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
				double ph = phi(s, w, Math.pow(z, rho) * Math.exp(x * sigma));
				p += ph * wp;
			}
			if (Math.abs(p - goal) < bestPhiError) {
				bestPhiError = Math.abs(p - goal);
				bestS = s;
			}
		}
		return bestS;
	}

	private double[] simulatedAnnealing(double[] k, double[] z, double[] guess) {

		Random rng = new Random(42);
		double[][] normal = new double[N][ROLLS];
		for (int n = 1; n < N; n++) {
			for (int r = 0; r < ROLLS; r++) {
				normal[n][r] = Math.pow(z[n], rho) * Math.exp(rng.nextGaussian() * sigma);
			}
		}

		double[] current = guess.clone();
		double[] wages = new double[N];
		double[][] rolls = new double[N][ROLLS];
		for (int n = 1; n < N; n++) {
			wages[n] = (1 - delta) * k[n] + z[n] * Math.pow(k[n], alpha);
			for (int r = 0; r < ROLLS; r++) {
				rolls[n][r] = phi(current[n], wages[n], normal[n][r]);
			}
		}
		double[] sums = new double[ROLLS];
		double score = 0.0;
		for (int r = 0; r < ROLLS; r++) {
			for (int n = 1; n < N; n++) {
				sums[r] += rolls[n][r] / N;
			}
			score += Math.abs(sums[r] - 1.0) / ROLLS;
		}
		// Save to restore
		double[][] saveRolls = new double[N][];
		for (int n = 0; n < N; n++) {
			saveRolls[n] = rolls[n].clone();
		}
		double[] saveSums = sums.clone();
		double saveScore = score;

		double bestScore = Double.MAX_VALUE;
		double[] bestResult = current.clone();
		System.out.println("Before optimization, average score of " + score + ".");

		for (int attempts = 0; attempts < NUMBER_OF_SA; attempts++) {
			if (attempts > 0) {
				score = saveScore;
				sums = saveSums.clone();
				for (int n = 1; n < N; n++) {
					rolls[n] = saveRolls[n].clone();
				}
				current = guess.clone();
			}
			long start = System.currentTimeMillis();
			long timeOut = start + TIME_PER_SA;
			double T0 = 0.1 * score / N;
			double[] newRolls = new double[ROLLS];
			while (System.currentTimeMillis() < timeOut) {
				double time = (System.currentTimeMillis() - start) * 1.0 / TIME_PER_SA;
				double temperature = T0 * Math.exp(-8 * time);
				int accepted = 0;
				for (int round = 0; round < 100; round++) {
					int i = 1 + rng.nextInt(N-1);
					double w = wages[i];
					double newS = current[i] + 0.02 * (1 - time) * (rng.nextBoolean() ? 1.0 : -1.0);
					if (newS < 0.10 || newS > 0.90) {
						continue;
					}
					double currentScore = 0.0;
					for (int r = 0; r < ROLLS; r++) {
						newRolls[r] = phi(newS, w, normal[i][r]);
						currentScore += Math.abs(sums[r] + (newRolls[r] - rolls[i][r]) / N - 1.0) / ROLLS;
					}
					if (rng.nextDouble() <= Math.exp(-(currentScore - score) / temperature)) {
						accepted++;
						current[i] = newS;
						for (int r = 0; r < ROLLS; r++) {
							sums[r] += (newRolls[r] - rolls[i][r]) / N;
							rolls[i][r] = newRolls[r];
						}
						score = currentScore;
						if (score < bestScore) {
							bestScore = score;
							bestResult = current.clone();
						}
					}
				}
				if (time < 0.2) {
					if (accepted > 80) {
						T0 = 0.66 * T0;
					}
					if (accepted < 60) {
						T0 = 1.33 * T0;
					}
				}
				if (time > 0.8) {
					if (accepted > 10) {
						T0 = 0.66 * T0;
					}
					if (accepted == 0) {
						T0 = 1.33 * T0;
					}
				}
			}
			System.out.println("End of attempt best score " + bestScore + ".");
		}
		System.out.println("After optimization, average score of " + bestScore + ".");
		double[] result = new double[N];
		for (int n = 1; n < N; n++) {
			result[n] = bestResult[n] * wages[n];
		}
		return result;
	}

	private double phi(double s, double w, double z) {
		
		double nc = (1 - delta) * w * (1 - s) + z * Math.pow(w * (1 - s), alpha);
		return beta * (1 - delta + alpha * z * Math.pow(w * (1 - s), alpha - 1)) * Math.pow(s * w / nc, eta);
	}

	public int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {

		System.out.println("beta = " + beta);
		System.out.println("eta = " + eta);
		System.out.println("alpha = " + alpha);
		System.out.println("delta = " + delta);
		System.out.println("rho = " + rho);
		System.out.println("sigma = " + sigma);
		this.t = 0;
		this.beta = beta;
		this.eta = eta;
		this.alpha = alpha;
		this.delta = delta;
		this.rho = rho;
		this.sigma = sigma;
		this.N = N;
		this.T = T;
		return 0;
	}

}
