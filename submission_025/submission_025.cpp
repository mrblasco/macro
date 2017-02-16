#include <vector>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iostream>
#include <time.h>
#include <emmintrin.h> 
#include <pmmintrin.h> 

#include <cassert>
#include <sys/time.h>


using namespace std;



#define VC		vector
#define VS		VC <string>
#define VI		VC <int>
#define VD		VC <double>
#define VVD		VC <VD>
#define S		size()
#define PB		push_back

	double getTime() {timeval tv;gettimeofday(&tv, NULL);return tv.tv_sec+tv.tv_usec*1e-6;}
	
	VS split(string s, char c) {
		VS st; int p=0,np; 
		while (np=s.find(c,p),np>=0) {
			//if (np!=p) 
			st.PB(s.substr(p,np-p));
			p=np+1;
		}
		if (p<s.S) st.PB(s.substr(p)); 
		return st;
	}

	//Mersenne twister
	struct RND {
	    int index;
	    unsigned int MT[624];
		
		RND(int seed=1) {init(seed);}
	    
	    void init(int seed=1) {
	        index=0;
	        MT[0]=seed;
	        for (int i=1;i<624;i++) MT[i]=(1812433253UL*(MT[i-1]^(MT[i-1]>>30))+i);
	    }
	    
	    void generate() {
	        const unsigned int MULT[2] = {0, 2567483615UL};
	        for (int i=0;i<227;i++) {
	            unsigned int y=(MT[i]&0x8000000UL)+(MT[i+1]&0x7FFFFFFFUL);
	            MT[i]=MT[i+397]^(y>>1);
	            MT[i]^=MULT[y&1];
	        }
	        for (int i=227;i<623;i++) {
	            unsigned int y=(MT[i]&0x8000000UL)+(MT[i+1]&0x7FFFFFFFUL);
	            MT[i]=MT[i-227]^(y>>1);
	            MT[i]^=MULT[y&1];
	        }
	        unsigned int y=(MT[623]&0x8000000UL)+(MT[0]&0x7FFFFFFFUL);
	        MT[623]=MT[623-227]^(y>>1);
	        MT[623]^=MULT[y&1];
	    }
	    
	    unsigned int rand() {
			if (index==0) generate();
	        
	        unsigned int y=MT[index];
	        y ^= y >> 11;
	        y ^= y <<  7 & 2636928640UL;
	        y ^= y << 15 & 4022730752UL;
	        y ^= y >> 18;
	        index=(index==623) ? 0 : index + 1;
	        return y;
	    }
	};
	
	struct RFparam {
		int RFeatures;
		int RPositions;
		int NodeSize;
		int Ntrees;
		double RFtime;
		int RFimportance;
		RFparam() {RFimportance=-1;}
	} RFp;

	struct Node {
		double value,result;
		int level,feature;
		int left,right;	
		Node() {value=result=0.;level=feature=left=right=-1;}
	};
	
	

//------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////

	//RandomForest RAF;
	VVD trainData,trainData0;
	VVD testData;
	VD results;
	double aver[5][15];int naver[5][15];

	int SI,SCE;
	bool frs;
	double tim0,endtim,totaltim,  sc;
	
	VS maxkod;
	
	struct fkods {
		string kod;
		int num;
	} fkod0;
	vector <vector <struct fkods> > fkod; 

//	vector < vector <int> > kod(30);vector <int> kod0;
//	vector < vector <int> > dat(30);vector <int> dat0;
//	int nuk[30];
	
	RND r(0);
		
#define NA -1000.
	
	struct Cmp1{bool operator()(fkods a, fkods b) { 
		if (a.num==b.num)return a.kod<b.kod;
		return a.num>b.num;
		} };





		int N,T;
		double beta,eta,alpha,delta,rho,sigma;	

		VVD AKt,AZt,ACt;
		int tt;
		double tottim,tim1,tim2;
	
//------------------------------------------------------------------------------	
	struct NationSave {

		int SetEconomyParameters(double beta1, double eta1, double alpha1, double delta1, double rho1, double sigma1, int N1, int T1){
			beta=beta1;eta=eta1;alpha=alpha1;delta=delta1;rho=rho1;sigma=sigma1;N=N1;T=T1;
			tt=0;
			static int si=0;si++;
			if (si==1)
			{tim0=getTime();tottim=0;}
		}

//------------------------------------------------------------------------------
		
		vector <double> ConsumptionDecisionRule(vector <double> Kt, vector <double> Zt){
			tim1=getTime();//if (tt==0)cerr<<"   tottim="<<tottim<<endl;
			if (tt==0){AKt.clear();AZt.clear();ACt.clear();}

			AKt.PB(Kt);AZt.PB(Zt);
			
			VD W(Kt.S);
			for (int i=0;i<Kt.S;i++)W[i]=(1-delta)*Kt[i]+Zt[i]*pow(Kt[i],alpha);
			
			VD Ct(Kt.S),CCt(Kt.S);
			
			double eert,meert=-1e9;
			if (tt==0){for (int i=0;i<Kt.S;i++)Ct[i]=W[i]/3;goto e0;}//10-967439
			if (tt==T){for (int i=0;i<Kt.S;i++)Ct[i]=W[i];goto e0;}

//if (tottim>-100){		
//			for (int i=0;i<Kt.S;i++){
//				double a=1-delta+Zt[i]*alpha*pow(Kt[i],alpha-1);
//				Ct[i]=ACt[tt-1][i]*pow(beta*a,1./eta);
//				Ct[i]=min(Ct[i],W[i]*.5);  
//			}
//			goto e0;
//		}

			for (int k=4;k<=9;k++){		//4-958	
				for (int i=0;i<Kt.S;i++){
					double a=1-delta+Zt[i]*alpha*pow(Kt[i],alpha-1);
					Ct[i]=ACt[tt-1][i]*pow(beta*a,1./eta);
					Ct[i]=min(Ct[i],W[i]*.1*k);  
				}
				eert=0;
				for (int i=0;i<Kt.S;i++){
						double eerti=beta*pow(ACt[tt-1][i]/Ct[i],eta)*(1-delta+Zt[i]*alpha*pow(Kt[i],alpha-1))-1;
						eert+=eerti;
				}
				if (eert>meert){meert=eert;for (int i=0;i<Kt.S;i++)CCt[i]=Ct[i];}
			}
			for (int i=0;i<Kt.S;i++)Ct[i]=CCt[i];

			e0:;
			
//			if (tt==0)for (int i=0;i<Kt.S;i++)Ct[i]=W[i]/3;
//			else if (tt<T){
//				for (int i=0;i<Kt.S;i++){
//					double a=1-delta+Zt[i]*alpha*pow(Kt[i],alpha-1);
//					Ct[i]=ACt[tt-1][i]/pow(beta*a,1./eta);
//				}					
//				for (int i=0;i<Kt.S;i++)Ct[i]=min(Ct[i],W[i]*.5);  
//			}
//			else 
//			for (int i=0;i<Kt.S;i++)Ct[i]=W[i];
			
			
			
			ACt.PB(Ct);
			
			for (int i=0;i<Kt.S;i++)if (Ct[i]>W[i])cerr<<tt<<"-"<<i<<": "<<Ct[i]<<"  "<<W[i]<<"     "<<endl;
			
			tt++;
			
			tim2=getTime();
			tottim+=tim2-tim1;
			return Ct;
			
		}
	
	
	};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////