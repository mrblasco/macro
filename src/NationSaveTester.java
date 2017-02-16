import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;
import java.io.File;
import java.io.FileOutputStream;
import java.io.DataOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


public class NationSaveTester
{
	static String mExecCommand = null;
	static String mTestFile = "./economies.csv";
	static String mShocksFile = "./shocks_system.csv";
	static String mOutputFile = "./competitor.csv";
	static int mTestOffset = 0;
	static int mNumTests = -1;
	static boolean mDebug = true;
	static Process mSolution = null;

	public void printMessage(String s)
	{
		if(mDebug)
		{
			System.out.println(s);
		}
	}

	public void doExec() throws Exception
	{

        // Send number of simulations N and number of time periods T
        final int N = 50;
        final int T = 2000;
        final double K0 = 1;
        final double Z0 = 1;

        // Read test cases file
        List<double[]> testCases = new ArrayList<double[]>();
        printMessage("Reading test case file: " + mTestFile + ".");
        BufferedReader br = new BufferedReader(new FileReader(mTestFile));
        br.readLine();
        while (true)
        {
            String s = br.readLine();
            if (s == null)
            {
                break;
            }
            String[] tokens = s.split(",");
            if(tokens.length != 6)
            {
            	printMessage("Incorrent number of tokens for test case " + s + ". Expected 6 got " + tokens.length + ".");
            	printMessage("Score = 0");
            	return;
            }

            double beta = Double.parseDouble(tokens[0]);
            double eta = Double.parseDouble(tokens[1]);
            double alpha = Double.parseDouble(tokens[2]);
            double delta = Double.parseDouble(tokens[3]);
            double rho = Double.parseDouble(tokens[4]);
            double sigma = Double.parseDouble(tokens[5]);

            testCases.add(new double[] { beta, eta, alpha, delta, rho, sigma });
        }

        // Read shocks
        double[][] shocks = new double[T][N];
        BufferedReader br_shocks = new BufferedReader(new FileReader(mShocksFile));
        for(int t = 0; t < T; t++)
        {
            String s = br_shocks.readLine();
            if (s == null) break;

            String[] tokens = s.split(",");
            if(tokens.length < N)
            {
            	printMessage("Incorrent number of shocks for time " + t + ". Expected " + N + " got " + tokens.length + ".");
            	printMessage("Score = 0");
            	return;
            }

            for(int n = 0; n < N; n++) shocks[t][n] = Double.parseDouble(tokens[n]);
        }

        Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mOutputFile), "utf-8"));
        writer.write("# Economy_id, simulation_id, period, C, Z, K, EulerResidual\n");
        Writer writer_log = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mOutputFile + ".log"), "utf-8"));
        writer_log.write("# Economy_id, score\n");

        // Execute test cases
        printMessage("Executing " + testCases.size() + " test cases.");

        printMessage("\t" + "Settings: "
            + "N: " + N + ", T: " + T
            + ", K0: " + K0 + ", Z0: " + Z0);

        double score = 0;
        int last_test = testCases.size();
        if (mNumTests >= 0 && mTestOffset + mNumTests < testCases.size()) last_test = mTestOffset + mNumTests;
        for (int test = mTestOffset; test < last_test; test++)
        {
        	printMessage("Test case " + (test + 1) + ".");

            NationSave predictor = new NationSave();

            boolean hasError = false;

        	// Send economy parameters
        	double beta = testCases.get(test)[0];
            double eta = testCases.get(test)[1];
            double alpha = testCases.get(test)[2];
            double delta = testCases.get(test)[3];
            double rho = testCases.get(test)[4];
            double sigma = testCases.get(test)[5];

            predictor.SetEconomyParameters(beta, eta, alpha, delta, rho, sigma, N, T);

			// Consumption
			double[][] C = new double[T + 1][N];
			// Capital
			double[][] K = new double[T + 1][N];
			// Wage
			double[][] W = new double[T + 1][N];
			// Shock
			double[][] Z = new double[T + 1][N];
			// Euler residuals
			double[][] eer = new double[T + 1][N];

			// Set initial capital
			Arrays.fill(K[0], K0);
			// Set initial shock
			Arrays.fill(Z[0], Z0);

			double testScore = 0.0;
			for(int t = 0; t <= T; t++)
			{
				for(int n = 0; n < N; n++)
				{
					// Update shock
					if(t != 0) //Already set for t = 0
					{
						double eps_t = sigma * shocks[t - 1][n];
						Z[t][n] = Math.exp(rho * Math.log(Z[t - 1][n]) + eps_t);
					}
					// Compute wage
					W[t][n] = (1 - delta) * K[t][n] + Z[t][n] * Math.pow(K[t][n], alpha);
				}
				if(t != T)
				{

					C[t] = predictor.ConsumptionDecisionRule(K[t], Z[t]);

					for(int i = 0; i < N; i++)
					{
						//Sanity check C[t][i]
						if(C[t][i] < 0)
						{
							printMessage("Consumption cannont be negative.");
							printMessage("Score = 0");
							testScore = -1;
							hasError = true;
							break;
						}
						if(C[t][i] > W[t][i])
						{
							printMessage("Consumption cannont exceed wage.");
							printMessage("Score = 0");
							testScore = -2;
							hasError = true;
							break;
						}
						if(Double.isNaN(C[t][i]))
						{
							printMessage("Consumption must be a real number.");
							printMessage("Score = 0");
							testScore = -3;
							hasError = true;
							break;
						}

						//Update capital for next period
						K[t + 1][i] = W[t][i] - C[t][i];
					}
				}
				else
				{
					//On last time period, consume all of the wage
					C[t] = W[t];
				}
/*
                if (t < 3) {
                    if (t == 0) {
                        printMessage(
                                "[" + "t" + "]"
                            + " \t" + "K[t][n]"
                            + " \t" + "Z[t][n]"
                            + " \t" + "W[t][n]"
                            + " \t" + "C[t][n]"
                            );
                    }
                    for (int n = 0; n < 4; n++) {
                        printMessage(
                                "[" + t + "]"
                            + " \t" + K[t][n]
                            + " \t" + Z[t][n]
                            + " \t" + W[t][n]
                            + " \t" + C[t][n]
                            );
                    }
                }
//*/
                if (hasError) break;

				if (t > 0)
                {
                    for(int i = 0; i < N; i++)
                    {
                        eer[t-1][i] += beta * Math.pow(C[t-1][i] / C[t][i], eta) * (1 - delta + Z[t][i] * alpha * Math.pow(K[t][i], alpha - 1)) - 1;
                    }
                }
			}

			// Scoring
            if (! hasError) {
                for (int t = 0; t < T; t++) {
                    double eer_t = 0;
                    for(int i = 0; i < N; i++)
                    {
                        eer_t += eer[t][i];
                    }
                    testScore += Math.abs(eer_t) / N;
                }
                testScore = 1000000 / (1 + testScore / T);
                score += testScore;
            }


			printMessage(Double.toString(testScore));

            if (! hasError) {
                printMessage("Saving...");
                // fout << "# Economy_id, simulation_id, period, C, Z, K, EulerResidual" << endl;
                for (int n = 0; n < N; n++) {
                    for (int t = 0; t <= T; t++) {
                        String str_out = "" + test + "," + n + "," + t + ","
                            + C[t][n] + "," + Z[t][n] + "," + K[t][n] + ",";
                        if (t < T) str_out += eer[t][n];
                        writer.write(str_out + "\n");
                    }
                }
            }

            writer_log.write("" + test + ",\t" + testScore + "\n");

        }
        score /= testCases.size();
        System.out.println("Score = " + score);
        writer.close();
        writer_log.close();
	}

	public static void main(String[] args)
	{
		for(int i = 0; i < args.length; i++)
		{
			if (args[i].equals("-eco"))
			{
                mTestFile = args[++i];
			}
			else if (args[i].equals("-shocks"))
			{
                mShocksFile = args[++i];
			}
			else if (args[i].equals("-out"))
			{
                mOutputFile = args[++i];
			}
			else if (args[i].equals("-test_offset"))
			{
                mTestOffset = Integer.parseInt(args[++i]);;
			}
			else if (args[i].equals("-num_tests"))
			{
                mNumTests = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-silent"))
			{
                mDebug = false;
			}
			else
			{
				System.out.println("WARNING: unknown argument " + args[i] + ".");
			}
		}

		try
		{
			new NationSaveTester().doExec();
		}
		catch (Exception e)
		{
			System.out.println("FAILURE: " + e.getMessage());
            e.printStackTrace();
            mSolution.destroy();
		}
	}

	class ErrorStreamRedirector extends Thread
	{
        public BufferedReader mReader;

        public ErrorStreamRedirector(InputStream is)
        {
        	mReader = new BufferedReader(new InputStreamReader(is));
        }

        public void run()
        {
            while (true)
            {
                String s;
                try
                {
                    s = mReader.readLine();
                }
                catch (Exception e)
                {
                    // e.printStackTrace();
                    return;
                }
                if (s == null) {
                    break;
                }
                System.out.println(s);
            }
        }
    }

}
