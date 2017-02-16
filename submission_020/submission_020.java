import java.lang.*;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class NationSave {
static double _beta;
static double _eta;
static double _alpha;
static double _delta;
static double _rho;
static double _sigma;
static int _N;
static int _T;
static int _tCount;
static boolean _firsttime = true;
static double[] _C;
static double[] _Cpre;
static int _N0;
static double _rmax;
static double _rmid;
static double _rmin;
static double _dmax;
static double _dmin;
static double _cmin;
static double _wmin;
static double _rmax0;
static double _rmid0;

static int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
	_beta = beta;
	_eta = eta;
	_alpha = alpha;
	_delta = delta;
	_rho = rho;
	_sigma = sigma;
	_N = N;
	_T = T;
	_tCount = 0;
	_rmax = 0.55770225;
	_rmid = 0.5;
	_rmin = 0.45;
	_dmax = 10;
	_dmin = 0.500001;
	_cmin = 0.0001;
	_wmin = 0.00001;

	_rmax0 = 0.55770225;
	_rmid0 = 0.5;
	//_rmax = adjust_r(_rmax0,1.5801275);
	//_rmax = adjust_r(_rmax0,1.855); // pow(1-_delta,1/_eta);
	_rmax = adjust_r(_rmax0,0.48222225); // sigmoid(_eta/(1-_delta));
	//_rmax = adjust_r(_rmax0,0.75); // sigmoid(1/(1-delta))
	//_rmax = adjust_r(_rmax0,2.255); // sigmoid(1-delta)
	//_rmid = adjust_r(_rmid0,2);

	if (_firsttime) {
		_firsttime = false;
		_C = new double[_N];
		_Cpre = new double[_N];
		_N0 = _N;
	} else {
		if (_N > _N0) {
			// gc will take care of the previously allocated memory??
			_C = null; System.gc(); // not necessary?? 
			_C = new double[_N];
			_Cpre = new double[_N];
			_N0 = _N;
		}
	}
	return 0;
}


static double sigmoid(double x) {
	double s = 1/(1+Math.pow(Math.E,-x));
	s -= 0.5;
	return s;
}

static double adjust_r(double r, double x) {
	//double x = 1.5801275; // for .100 and no eta
	//double x = -12.20975;
	//double scaler = Math.pow(1-_delta,1/_eta);
	//double scaler = Math.pow(1-_delta,1/_eta);
	double scaler = _eta/(1-_delta);
	double s = sigmoid(x*scaler);
	//double s = Math.pow(scaler,x);
	double dy = Math.min(1-r,r);
	r += dy * s;
	return r;
}


static double bound_check(double c, double w) {
	if (c > _rmax*w) {
		c = _rmax * w;
	} else if (w<=_wmin) {
		c = _rmid*w;
	}  
    if(Double.isNaN(c)) {
    	c = 0;
    }
    if (c>=w) {
     	//c = w;
    	c = _rmax * w; // c=w => k[t+1]=0 and may cause scoring problems score=Nan
    }
    //if (c<=0.1) {c = 0;} 

	return c;
}

static double calc_decay_forward(double k, double z) {
    double v = k>0 ? _beta*(1-_delta+z*_alpha*Math.pow(k,_alpha-1)) : _beta*(1-_delta);
    double decay = v>0 ? Math.pow(v,1/_eta) : 0;
    return decay;
}

static double optimal_consumption_perT(double w, double cpre, double z, double k) {
    //double c = _cmin;
    double c = _rmid*w;
    if (w > _wmin) {
    	double decay = calc_decay_forward(k,z);
    	//c = cpre * decay;
    	if (decay <_dmin) {
    		//c = _rmax*w;
    		c = cpre>0 ? cpre*_dmin : w*_rmin;
    	} else if (decay>_dmax) {
    		//c = _rmax*w;
    		c = cpre>0 ? cpre*_dmax : w*_rmax;
    	} else {
	    	c = cpre>0 ? cpre*decay : c;
    	}
    	c = bound_check(c,w);
		//System.err.println(String.format("optimal_consumption_perT %f %f %f %f => c=%f r=%f",w,cpre,z,k,c,w>0?c/w:0));
	    //System.err.println(String.format("optimal_consumption_perT %.3f",c));
	}
   return c;
}


static double[] ConsumptionDecisionRule(double[] K, double[] Z) {
	int t_remain = _T - _tCount;
	double r = _rmax;
	//if (_tCount==0 || t_remain <= 4) {
	//if (_tCount==0 || t_remain <= 1) {
	if (_tCount==0) {
	//if (false) {
		if (t_remain == 4) {
			//r = 0.53;
			r = 0.515;
		} else if (t_remain == 3) {
			//r = 0.58;
			r = 0.54;
		} else if (t_remain == 2) {
			//r = 0.65;
			r = 0.58;
		} else if (t_remain == 1) {
			//r = 0.75;
			r = 0.68;
		} 
		for (int n=0; n<_N; n++) {
			double w = K[n]<=0 ? 0 : (1 - _delta) * K[n] + Z[n] * Math.pow(K[n], _alpha);
			//if (w<0) { w=0; }
			double c = w>_wmin ? r*w : 0;
			//c = bound_check(c,w);
			_C[n] = c;
			_Cpre[n] = _C[n];
			//System.err.println(String.format("ConsumptionDecisionRule %d z=%f k=%f => w=%f c=%f r=%f",n,Z[n],K[n],w,_C[n],w>0 ? _C[n]/w : 0));
		}
	} else {
		for (int n=0; n<_N; n++) {
			double w = K[n] > 0 ? (1 - _delta) * K[n] + Z[n] * Math.pow(K[n], _alpha) : 0;
			//_C[n] =  r*w;
			double c = optimal_consumption_perT(w,_Cpre[n],Z[n],K[n]);
			_C[n] = c;
			_Cpre[n] = _C[n];
			//_C[n] = r*w;
			//System.err.println(String.format("ConsumptionDecisionRule %d z=%f k=%f => w=%f c=%f r=%f",n,Z[n],K[n],w,_C[n],w>0 ? _C[n]/w : 0));
		}
	}
	_tCount++;
	return _C;
}


public static void main(String[] args) {

	Scanner sc = new Scanner(System.in);
	int X = sc.nextInt();
	for (int x=0; x<X; x++) {
		double beta = sc.nextDouble();
		double eta = sc.nextDouble();
		double alpha = sc.nextDouble();
		double delta = sc.nextDouble();
		double rho = sc.nextDouble();
		double sigma = sc.nextDouble();
		int N = sc.nextInt();
		int T = sc.nextInt();

		SetEconomyParameters(beta, eta, alpha, delta, rho, sigma, N, T);

	 	double [] K = new double[N];
	 	double [] Z = new double[N];
	 	// don't need to feed data to tests for t=T
		for (int t=0; t<T; t++) {
			for (int n=0; n<N; n++) {
				double k  = sc.nextDouble();
				K[n] = k;
			}
			for (int n=0; n<N; n++) {
				double z = sc.nextDouble();
				Z[n] = z;
			}
			//System.err.println(String.format("t=%d/%d",t,T));
			ConsumptionDecisionRule(K,Z);
			for (int n=0; n<N; n++) {
				System.out.println(_C[n]);
			}
			System.out.flush();
		}
	}

}
}


