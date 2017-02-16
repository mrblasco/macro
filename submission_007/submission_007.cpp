#define NDEBUG

#define GV_JS

#undef __STRICT_ANSI__

#include <string>
#include <cstdlib>
#include <cstdio>
#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <functional>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cassert>

using namespace std;

#ifdef NDEBUG
#undef assert
#define assert(e)
#else

class MemoStacktrace {
	static char buf[65536];
	static vector<int> stack;
public:
	MemoStacktrace(const char * func, const char * file, int line) {
		int now = stack.empty() ? 0 : stack.back();
		int now2 = now + sprintf(buf + now, " ... %s: %s, %d\n", func, file, line);
		if(32768<now2) {
			now2 = now;
			fprintf(stderr, "ERROR in MemoStacktrace ... stacktrace is overflow!\n");
			fflush(stderr);
		}
#ifdef STACKTRACE_FORCE_PRINT_FOR_DEBUG
		char sp[256];
		for(int i=0; i<max(40, (int)stack.size()); ++i) {
			sp[i] = ' ';
		}
		sp[stack.size()] = '\0';
		fprintf(stderr, "!%s%s: %s, %d {\n", sp, func, file, line);
		fflush(stderr);
#endif
		stack.push_back(now2);
	}
	~MemoStacktrace() {
		if(stack.empty()) {
			fprintf(stderr, "ERROR in MemoStacktrace ... why stack is empty!?\n");
			fflush(stderr);
		}
		else {
			stack.pop_back();
		}
#ifdef STACKTRACE_FORCE_PRINT_FOR_DEBUG
		char sp[256];
		for(int i=0; i<max(40, (int)stack.size()); ++i) {
			sp[i] = ' ';
		}
		sp[stack.size()] = '\0';
		fprintf(stderr, "!%s}\n", sp);
		fflush(stderr);
#endif
	}
	static void print() {
		int now = stack.empty() ? 0 : stack.back();
		if(now) {
			buf[now] = '\0';
			fprintf(stderr, "[StackTrace]\n");
			fwrite(buf, 1, now, stderr);
			fflush(stderr);
		}
	}
};
char MemoStacktrace::buf[65536];
vector<int> MemoStacktrace::stack;
#if defined(STACKTRACE_FOR_DEBUG) || defined(STACKTRACE_FORCE_PRINT_FOR_DEBUG)
#else
#endif

inline void myAbort() {
	MemoStacktrace::print();
	abort();
}
inline void myAssert(bool v, const char * func, const char * file, int line, const char * expr) {
	if(__builtin_expect(!v, 0)) {
		fprintf(stderr, "Assertion failed: (%s) in (%s, %d: %s).\n", expr, file, line, func);
		fflush(stderr);
		myAbort();
	}
}
#undef assert
#define assert(e) myAssert(e, __PRETTY_FUNCTION__, __FILE__, __LINE__, #e)
#endif

long long getTestcaseSeed(long long set=0) {
	static long long seed = 0;
	if(set) {
		seed = set;
	}
	return seed;
}

template<typename T, typename TyCompare = less<T> > class RankingVector {
	vector<T> vec;
	vector<int> ranking;
	vector<int> revRanking;
	inline void pushRanking(int idx) {
		assert(vec.size()==revRanking.size());
		assert(idx<(int)vec.size());
		const T & value = vec[idx];
		int now = ranking.size();
		ranking.emplace_back();
		if(2<=now) {
			TyCompare cmp;
			int clz = __builtin_clz(now+2);

			if(now<((3<<(30-clz))-2)) {
				int nxt = (now-2)>>1;
				if(cmp(vec[ranking[nxt]], value)) {
					do {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					} while(cmp(vec[ranking[nxt]], value));
				}
				else {
					int nxt = (now+(1<<(30-clz))-2)>>1;
					while(cmp(value, vec[ranking[nxt]])) {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					}
				}
			}
			else {
				int nxt = (now-2)>>1;
				if(cmp(value, vec[ranking[nxt]])) {
					do {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					} while(cmp(value, vec[ranking[nxt]]));
				}
				else {
					int nxt = now-(1<<(30-clz));
					while(cmp(vec[ranking[nxt]], value)) {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					}
				}
			}
		}
		else if(now==1) {
			TyCompare cmp;
			if(cmp(vec[ranking[0]], value)) {
				revRanking[ranking[1] = ranking[0]] = 1;
				now = 0;
			}
		}
		revRanking[ranking[now] = idx] = now;
	}
	inline void adjustDownRanking(int now, int idx, const T & value) {
		TyCompare cmp;
		while(true) {
			int nxt1 = now + now + 2;
			int nxt2 = nxt1 + 1;
			if(nxt2<(int)ranking.size()) {
				if(cmp(vec[ranking[nxt1]], vec[ranking[nxt2]])) {
					if(cmp(value, vec[ranking[nxt2]])) {
						revRanking[ranking[now] = ranking[nxt2]] = now;
						now = nxt2;
					}
					else {
						revRanking[ranking[now] = idx] = now;
						return;
					}
				}
				else {
					if(cmp(value, vec[ranking[nxt1]])) {
						revRanking[ranking[now] = ranking[nxt1]] = now;
						now = nxt1;
					}
					else {
						revRanking[ranking[now] = idx] = now;
						return;
					}
				}
			}
			else if(nxt1<(int)ranking.size()) {
				if(cmp(value, vec[ranking[nxt1]])) {
					revRanking[ranking[now] = ranking[nxt1]] = now;
					now = nxt1;
				}
				else {
					revRanking[ranking[now] = idx] = now;
					return;
				}
			}
			else {
				break;
			}
		}
		int clz = __builtin_clz(now+2);

		int nxt = now+(1<<(30-clz));
		if((int)ranking.size()<=nxt) {
			nxt = (nxt-2)>>1;
		}
		if(0<=nxt) {
			while(cmp(value, vec[ranking[nxt]])) {
				revRanking[ranking[now] = ranking[nxt]] = now;
				now = nxt;
				if(now<2) {
					break;
				}
				nxt = (now-2)>>1;
			}
		}
		revRanking[ranking[now] = idx] = now;
	}
	inline void adjustUpRanking(int now, int idx, const T & value) {
		TyCompare cmp;
		while(true) {
			int nxt1 = now + now + 2;
			int nxt2 = nxt1 + 1;
			if(nxt2<(int)ranking.size()) {
				if(cmp(vec[ranking[nxt2]], vec[ranking[nxt1]])) {
					if(cmp(vec[ranking[nxt2]], value)) {
						revRanking[ranking[now] = ranking[nxt2]] = now;
						now = nxt2;
					}
					else {
						revRanking[ranking[now] = idx] = now;
						return;
					}
				}
				else {
					if(cmp(vec[ranking[nxt1]], value)) {
						revRanking[ranking[now] = ranking[nxt1]] = now;
						now = nxt1;
					}
					else {
						revRanking[ranking[now] = idx] = now;
						return;
					}
				}
			}
			else if(nxt1<(int)ranking.size()) {
				if(cmp(vec[ranking[nxt1]], value)) {
					revRanking[ranking[now] = ranking[nxt1]] = now;
					now = nxt1;
				}
				else {
					revRanking[ranking[now] = idx] = now;
					return;
				}
			}
			else {
				break;
			}
		}
		int clz = __builtin_clz(now+2);

		int nxt = now+now-(1<<(31-clz))+2;
		if((int)ranking.size()<=nxt) {
			nxt = (nxt-2)>>1;
		}
		if(0<=nxt) {
			while(cmp(vec[ranking[nxt]], value)) {
				revRanking[ranking[now] = ranking[nxt]] = now;
				now = nxt;
				if(now<2) {
					break;
				}
				nxt = (now-2)>>1;
			}
		}
		revRanking[ranking[now] = idx] = now;
	}
	inline void adjustRanking(int idx) {
		assert(vec.size()==revRanking.size());
		assert(0<=idx && idx<(int)vec.size());
		const T & value = vec[idx];
		int now = revRanking[idx];
		int clz0 = __builtin_clz(now+2);

		TyCompare cmp;
		if(now<((3<<(30-clz0))-2)) {
			if(2<=now) {
				int nxt = (now-2)>>1;
				if(cmp(vec[ranking[nxt]], value)) {
					do {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					} while(cmp(vec[ranking[nxt]], value));
					revRanking[ranking[now] = idx] = now;
					return;
				}
			}
			adjustDownRanking(now, idx, value);
		}
		else {
			if(2<=now) {
				int nxt = (now-2)>>1;
				if(cmp(value, vec[ranking[nxt]])) {
					do {
						revRanking[ranking[now] = ranking[nxt]] = now;
						now = nxt;
						if(now<2) {
							break;
						}
						nxt = (now-2)>>1;
					} while(cmp(value, vec[ranking[nxt]]));
					revRanking[ranking[now] = idx] = now;
					return;
				}
			}
			adjustUpRanking(now, idx, value);
		}
	}
public:
	inline void clear() {
		vec.clear();
		ranking.clear();
		revRanking.clear();
	}
	inline void push_back(const T & o) {
		assert(vec.size()==revRanking.size());
		assert(vec.size()==ranking.size());
		int idx = vec.size();
		vec.push_back(o);
		revRanking.emplace_back();
		pushRanking(idx);
	}
	inline const T & get(int idx) {
		return vec[idx];
	}
	inline void set(int idx, const T & o) {
		vec[idx] = o;
		adjustRanking(idx);
	}
	inline int top() const {
		assert(0<ranking.size());
		return ranking[0];
	}
	inline int bottom() const {
		assert(0<ranking.size());
		return ranking.size()==1 ? ranking[0] : ranking[1];
	}
	inline bool empty() const {
		return vec.empty();
	}
	inline int size() const {
		return vec.size();
	}
};

#include <cstdarg>

const char * format(const char * fmt, ...) {
	static char buf[65536];
	va_list arg;
	va_start(arg, fmt);
	vsprintf(buf, fmt, arg);
	va_end(arg);
	return buf;
}

void info(const char * format, ...) {
	va_list arg;
	va_start(arg, format);
	vfprintf(stderr, format, arg);
	va_end(arg);
}

#define dump(...)

typedef unsigned int uint;
typedef unsigned int uint32;

#define arraysizeof(ARRAY) ((int)(sizeof(ARRAY) / sizeof(*(ARRAY))))

void gvUsing(function<void()> func) {
	func();
}

#include <sys/time.h>

class XsRandom {
	unsigned long long a;
	unsigned long long b;
public:
	inline XsRandom() : a(0x8a5cd789635d2dffULL), b(0x121fd2155c472f96ULL) {
	}
	inline XsRandom(const XsRandom & o) : a(o.a), b(o.b) {
	}
	inline unsigned long long next64() {
		unsigned long long c = a ^ (a<<23);
		a = b;
		b = c ^ a ^ (c>>18) ^ (a>>5);
		return b + a;
	}
	inline XsRandom(unsigned int seed) : a(0x8a5cd789635d2dffULL), b(0x121fd2155c472f96ULL) {
		seed = seed * 1234567891 + 521288629;
		unsigned long long a2 = a;
		unsigned long long b2 = b;
		while(seed) {
			next64();
			if(seed & 1) {
				a2 ^= a;
				b2 ^= b;
			}
			seed >>= 1;
		}
		a = a2;
		b = b2;
	}
	inline unsigned int next() {
		return (unsigned int)next64();
	}
	inline unsigned int next2() {
		unsigned int a = next();
		return (((unsigned long long)a + 1) * a)>>32;
	}
	inline int nextInt(int r) {
		assert(1<=r);
		return ((unsigned long long)next() * r)>>32;
	}
	inline double nextFloat() {
		return (next()+0.5) * (1.0/4294967296.0);
	}
	inline double next2Float() {
		return (next2()+0.5) * (1.0/4294967296.0);
	}
};

class PcgRandom {
	unsigned long long state;
public:
	inline PcgRandom() : state(0x8a5cd789635d2dffULL) {
	}
	inline PcgRandom(const PcgRandom & o) : state(o.state) {
	}
	inline unsigned int next() {
		unsigned long long old = state;
		state = old * 6364136223846793005ULL + 1442695040888963407ULL;
		unsigned int ret = ((old>>18)^old) >> 27;
		unsigned int rot = (old>>59);
		return (ret>>rot) | (ret<<-rot);
	}
	inline PcgRandom(unsigned int seed) : state(0x8a5cd789635d2dffULL) {
		seed = seed * 1234567891 + 521288629;
		unsigned long long st = state;
		while(seed) {
			next();
			if(seed & 1) {
				st ^= state;
			}
			seed >>= 1;
		}
		state = st;
	}
	inline unsigned int next2() {
		unsigned int a = next();
		return (((unsigned long long)a + 1) * a)>>32;
	}
	inline int nextInt(int r) {
		assert(1<=r);
		return ((unsigned long long)next() * r)>>32;
	}
	inline double nextFloat() {
		return (next()+0.5) * (1.0/4294967296.0);
	}
	inline double next2Float() {
		return (next2()+0.5) * (1.0/4294967296.0);
	}
};

class MTRandom {
	unsigned int state[624];
	int counter;
public:
	inline MTRandom(const MTRandom & o) {
		memcpy(this, &o, sizeof(MTRandom));
	}
	inline unsigned int next() {
		if(arraysizeof(state)<=counter) {
			unsigned int vv[] = { 0, 0x9908B0DF };
			int i;
			for(i=0; i<arraysizeof(state)-397; ++i) {
				unsigned int v = (state[i]&0x80000000) | (state[i+1]&0x7FFFFFFF);
				state[i] = state[i+397] ^ (v>>1) ^ vv[v&1];
			}
			for(; i<arraysizeof(state)-1; ++i) {
				unsigned int v = (state[i]&0x80000000) | (state[i+1]&0x7FFFFFFF);
				state[i] = state[i+397-arraysizeof(state)] ^ (v>>1) ^ vv[v&1];
			}
			{
				unsigned int v = (state[arraysizeof(state)-1]&0x80000000) | (state[0]&0x7FFFFFFF);
				state[arraysizeof(state)-1] = state[397] ^ (v>>1) ^ vv[v&1];
			}
			counter = 0;
		}
		unsigned int ret = state[counter];
		++counter;
		ret ^= ret>>11;
		ret ^= (ret<<7) & 0x9D2C5680;
		ret ^= (ret<<15) & 0xEFC60000;
		ret ^= ret>>18;
		return ret;
	}
	inline MTRandom(unsigned int seed = 4357) {
		for(int i=0; i<arraysizeof(state); ++i) {
			unsigned int v = seed & 0xFFFF0000;
			seed = 69069 * seed + 1;
			state[i] = v | ((seed>>16) & 0x0000FFFF);
			seed = 69069 * seed + 1;
		}
		counter = arraysizeof(state);
	}
	inline unsigned int next2() {
		unsigned int a = next();
		return (((unsigned long long)a + 1) * a)>>32;
	}
	inline int nextInt(int r) {
		assert(1<=r);
		return ((unsigned long long)next() * r)>>32;
	}
	inline double nextFloat() {
		return (next()+0.5) * (1.0/4294967296.0);
	}
	inline double next2Float() {
		return (next2()+0.5) * (1.0/4294967296.0);
	}
};

#if defined(USE_RANDOM_MT)
typedef MTRandom MyRandom;
#elif defined(USE_RANDOM_PCG)
typedef PcgRandom MyRandom;
#elif defined(USE_RANDOM_XS)
typedef XsRandom MyRandom;
#else
typedef XsRandom MyRandom;
#endif

MyRandom g_myRand;

static inline unsigned int myRand() {
	return g_myRand.next();
}
static inline void myRandInit(unsigned int seed = 0) {
	g_myRand = MyRandom(seed);
}

static inline int myRandInt(int r) {
	return g_myRand.nextInt(r);
}
static inline double myRand2Float() {
	return g_myRand.next2Float();
}

static inline double getTime2(){
	timeval tv;
	gettimeofday(&tv, 0);
	double result = tv.tv_sec + tv.tv_usec * 1e-6;
	return result;
}
static double g_startTime;
static double g_suspendTime = 0;
static inline void initTimeWithMyRandInit(){
	timeval tv;
	gettimeofday(&tv, 0);
	double result = tv.tv_sec + tv.tv_usec * 1e-6;
	myRandInit(tv.tv_usec);
	g_startTime = result;
	g_suspendTime = 0;
}
const double g_timeupSecBase = 9.8;
double g_timeupSec = g_timeupSecBase;

double watch(const char * format, function<void()> func, ...) {
	double start = getTime2();
	func();
	double spend = getTime2() - start;
	{
		va_list arg;
		va_start(arg, func);
		char buf[4096];
		int ret = vsnprintf(buf, sizeof(buf), format, arg);
		va_end(arg);
		if((int)sizeof(buf)<=ret) {
			buf[sizeof(buf)-1] = '\0';
		}
		info("%s: %fsec\n", buf, spend);
	}

	return spend;
}

double watch(function<void()> func) {
	return watch("time", func);
}

#include <numeric>

template<bool bigIsBetter>double searchConvex(function<double(double)> func, double a, double d, int n=20) {
	double b = (d-a)*0.381966+a;
	double c = (a-d)*0.381966+d;
	double bV = func(b);
	double cV = func(c);
	for(int i=0; i<n; ++i) {
		if(bigIsBetter ? cV<bV : bV<cV) {
			d = c;
			c = b;
			b = (d-a)*0.381966+a;
			cV = bV;
			bV = func(b);
		}
		else {
			a = b;
			b = c;
			c = (a-d)*0.381966+d;
			bV = cV;
			cV = func(c);
		}
	}
	if(bigIsBetter ? cV<bV : bV<cV) {
		return b;
	}
	else {
		return c;
	}
}

double paramA = 1.0;

double paramJ = 0.25;

double paramL = 1.0;

double paramN = 1.0;
double paramQ = 1.0;
double paramR = 8.0;
double paramU = 0.25;

void initParam(int kind = -1) {
	paramA = 0.9207642955;
	paramJ = 0.0822250266;
	paramL = 0.5670116868;
	paramN = 0.9920850012;
	paramQ = 1.3125476456;
	paramR = 8.2144599632;
	paramU = 0.1778701602;
}

static const int paramK = 1;

int g_count = 0;

double g_diffBaseA;
double g_diffBaseB;
double g_diffBaseAsmall;
double g_diffBaseBsmall;
double g_diffBaseAbig;
double g_diffBaseBbig;
double g_nextAB;

double funcFast(double v) {
  double diffA = g_diffBaseA - log(v);
  double diffB = g_diffBaseB - log(g_nextAB-v);
  return diffA*diffA+diffB*diffB;
}

class NationSave {
	double beta;
	double eta;
	double alpha;
	double delta;
	double rho;
	double sigma;
	int N;
	int T;
	int turn;
	vector<double> B;
	vector<double> targetRateList;
public:
	int SetEconomyParameters(double beta, double eta, double alpha, double delta, double rho, double sigma, int N, int T) {
		initTimeWithMyRandInit();
		initParam();
		++g_count;
		fprintf(stderr, "g_count: %d\n", g_count);
		dump(beta, eta, alpha, delta);
		dump(rho, sigma, N, T);
		this->beta = beta;
		this->eta = eta;
		this->alpha = alpha;
		this->delta = delta;
		this->rho = rho;
		this->sigma = sigma;
		this->N = N;
		this->T = T;
		turn = 0;
		{
			double lastRate = pow(beta*(1-delta*(1-alpha*pow(0.2, alpha-1))), 1/eta);
			static const double paramI = 0.00001;
			double paramJ2 = paramJ/lastRate;
			double v = 1.0;
			{

				targetRateList.push_back(v);
			}
			int t;
			for(t=1; t<paramK; ++t) {
				v = max(v*paramJ2, 0.001);
				targetRateList.push_back(v);
			}
			for(; t<T; ++t) {
				v = max(v*paramI, 0.001);
				targetRateList.push_back(v);
			}
			reverse(targetRateList.begin(), targetRateList.end());
			for(int t=0; t<T; ++t) {
				double sum = 0;
				for(int t2=t; t2<T; ++t2) {
					sum += targetRateList[t2];
				}
				sum += paramL;
				targetRateList[t] /= sum;
			}
		}
		return 0;
	}
	vector<double> ConsumptionDecisionRule(const vector<double> & K, const vector<double> & Z) {

		vector<double> C;
		if(turn==0) {
			double targetRate = targetRateList[turn];
			for(int i=0; i<(int)Z.size(); ++i) {
				double W = (1-delta) * K[i] + Z[i] * pow(K[i], alpha);
				C.push_back(W * targetRate);
			}
		}
		else if(turn+1==T) {
			vector<double> W;
			vector<double> now;
			vector<double> minNow;
			vector<double> Z2;
			vector<double> Z2small;
			vector<double> Z2big;
			vector<double> diffList;
			vector<double> errList;
			vector<double> err2List;
			vector<double> err3List;
			RankingVector<double> pena2;
			RankingVector<double> pena3;
			double diffSum = 0;
			double errSum = 0;
			double err2Sum = 0;
			double err3Sum = 0;
			double eSigma = exp(sigma*paramN);
			for(int i=0; i<(int)Z.size(); ++i) {
				double w = (1-delta) * K[i] + Z[i] * pow(K[i], alpha);
				W.push_back(w);
				C.push_back(B[i] * pow(beta*(1-delta+alpha*Z[i]*pow(K[i], alpha-1.0)), 1.0/eta));
				now.push_back(1.0);
				minNow.push_back(pow(C.back() / w, eta));
				Z2.push_back(pow(Z[i], rho));
				Z2small.push_back(Z2.back()/eSigma);
				Z2big.push_back(Z2.back()*eSigma);

				double aC = C.back();
				double aK = W.back() - aC;
				double aC2base = aC * pow(beta*(1-delta+alpha*Z2.back()*pow(aK, alpha-1.0)), 1.0/eta);
				double aW2 = (1-delta) * aK + Z2.back() * pow(aK, alpha);
				double aW2small = (1-delta) * aK + Z2small.back() * pow(aK, alpha);
				double aW2big = (1-delta) * aK + Z2big.back() * pow(aK, alpha);
				double eerA = beta * pow(aC / aW2, eta) * (1 - delta + Z2.back() * alpha * pow(aW2, alpha - 1)) - 1;
				double eerAsmall = beta * pow(aC / aW2small, eta) * (1 - delta + Z2small.back() * alpha * pow(aW2small, alpha - 1)) - 1;
				double eerAbig = beta * pow(aC / aW2big, eta) * (1 - delta + Z2big.back() * alpha * pow(aW2big, alpha - 1)) - 1;
				diffList.push_back(eerA);
				diffSum += eerA;
				double errAsmall = eerAsmall - eerA;
				double errAbig = eerAbig - eerA;
				double errA = errAsmall * errAsmall + errAbig * errAbig;
				double err2A = errAbig;
				double err3A = errAsmall;

				errList.push_back(errA);
				errSum += errA;
				err2List.push_back(err2A);
				err2Sum += err2A;
				err3List.push_back(err3A);
				err3Sum += err3A;
				pena2.push_back(eerA);
				pena3.push_back(minNow.back());
			}
			double score = diffSum*diffSum*paramQ + errSum*paramR + (diffSum+err2Sum+err3Sum)*(diffSum+err2Sum+err3Sum)*paramU;
			for(int i=0; i<100000; ++i) {

				int a, b;
				if(i&3) {
					a = myRandInt(now.size());
					b = myRandInt(now.size()-1);
					if(a<=b) {
						++b;
					}
				}
				else   {
					a = pena3.bottom();
					b = pena3.top();
				}

				double nowA = now[a];
				double nowB = now[b];
				uint c = myRand();
				if(true) {
					if(minNow[a]<nowA && (c&2)) {
						double d = myRand2Float() * (nowA-minNow[a]) * paramA;
						nowA -= d;
						nowB += d;
					}
					else {
						if(nowB<minNow[b]) {
							continue;
						}
						double d = myRand2Float() * (nowB-minNow[b]) * paramA;
						nowB -= d;
						nowA += d;
					}
				}
				double aC = C[a] * pow(nowA, -1.0/eta);
				double bC = C[b] * pow(nowB, -1.0/eta);
				double aK = W[a] - aC;
				double bK = W[b] - bC;

				double aW2 = (1-delta) * aK + Z2[a] * pow(aK, alpha);
				double aW2small = (1-delta) * aK + Z2small[a] * pow(aK, alpha);
				double aW2big = (1-delta) * aK + Z2big[a] * pow(aK, alpha);
				double bW2 = (1-delta) * bK + Z2[b] * pow(bK, alpha);
				double bW2small = (1-delta) * bK + Z2small[b] * pow(bK, alpha);
				double bW2big = (1-delta) * bK + Z2big[b] * pow(bK, alpha);
				if(aC<0.00001 || bC<0.00001 || aW2<0.00001 || bW2<0.00001) {
					continue;
				}

				double eerA = beta * pow(aC / aW2, eta) * (1 - delta + Z2[a] * alpha * pow(aK, alpha - 1)) - 1;
				double eerAsmall = beta * pow(aC / aW2small, eta) * (1 - delta + Z2small[a] * alpha * pow(aK, alpha - 1)) - 1;
				double eerAbig = beta * pow(aC / aW2big, eta) * (1 - delta + Z2big[a] * alpha * pow(aK, alpha - 1)) - 1;

				double eerB = beta * pow(bC / bW2, eta) * (1 - delta + Z2[b] * alpha * pow(bK, alpha - 1)) - 1;
				double eerBsmall = beta * pow(bC / bW2small, eta) * (1 - delta + Z2small[b] * alpha * pow(bK, alpha - 1)) - 1;
				double eerBbig = beta * pow(bC / bW2big, eta) * (1 - delta + Z2big[b] * alpha * pow(bK, alpha - 1)) - 1;

				double newDiffSum = diffSum - diffList[a] - diffList[b] + eerA + eerB;

				double errAsmall = eerAsmall - eerA;
				double errAbig = eerAbig - eerA;
				double errA = errAsmall * errAsmall + errAbig * errAbig;
				double err2A = errAbig;
				double err3A = errAsmall;

				double errBsmall = eerBsmall - eerB;
				double errBbig = eerBbig - eerB;
				double errB = errBsmall * errBsmall + errBbig * errBbig;
				double err2B = errBbig;
				double err3B = errBsmall;

				double newErrSum = errSum - errList[a] - errList[b] + errA + errB;
				double newErr2Sum = err2Sum - err2List[a] - err2List[b] + err2A + err2B;
				double newErr3Sum = err3Sum - err3List[a] - err3List[b] + err3A + err3B;

				double newSc = newDiffSum*newDiffSum*paramQ + newErrSum*paramR + (newDiffSum+newErr2Sum+newErr3Sum)*(newDiffSum+newErr2Sum+newErr3Sum)*paramU;
				if(newSc<score) {
					score = newSc;
					diffSum = newDiffSum;
					errSum = newErrSum;
					err2Sum = newErr2Sum;
					err3Sum = newErr3Sum;
					now[a] = nowA;
					now[b] = nowB;
					diffList[a] = eerA;
					diffList[b] = eerB;
					errList[a] = errA;
					errList[b] = errB;
					err2List[a] = err2A;
					err2List[b] = err2B;
					err3List[a] = err3A;
					err3List[b] = err3B;
					pena2.set(a, eerA);
					pena2.set(b, eerB);
					pena3.set(a, minNow[a]-nowA);
					pena3.set(b, minNow[b]-nowB);
				}
			}
			for(int i=0; i<(int)C.size(); ++i) {
				C[i] *= pow(now[i], -1.0/eta);
			}
			dump(turn, C[0]);
			dump(diffSum, errSum, err2Sum, err3Sum, score);
			for(int i=0; i<(int)C.size(); ++i) {
				if(W[i]<=C[i]) {
					C[i] = (W[i] / (1+T-turn));
				}
			}
		}
		else {
			vector<double> W;
			vector<double> now;
			vector<double> next;
			vector<double> minNow;
			vector<double> Z2;
			vector<double> pena;
			RankingVector<double> pena2;
			RankingVector<double> pena3;
			double targetRate = targetRateList[turn];
			for(int i=0; i<(int)Z.size(); ++i) {
				double w = (1-delta) * K[i] + Z[i] * pow(K[i], alpha);
				W.push_back(w);
				C.push_back(B[i] * pow(beta*(1-delta+alpha*Z[i]*pow(K[i], alpha-1.0)), 1.0/eta));
				now.push_back(1.0);
				next.push_back(1.0);
				minNow.push_back(pow(C.back() / w, eta));
				Z2.push_back(pow(Z[i], rho));

				double aC = C.back();
				double aK = W.back() - aC;
				double aC2base = aC * pow(beta*(1-delta+alpha*Z2.back()*pow(aK, alpha-1.0)), 1.0/eta);
				double aW2 = (1-delta) * aK + Z2.back() * pow(aK, alpha);
				double diffBaseA = (log(aW2 * targetRate) - log(aC2base))*eta;
				pena.push_back(diffBaseA*diffBaseA);
				pena2.push_back(diffBaseA);
				pena3.push_back(minNow.back());
			}
			for(int i=0; i<10000; ++i) {
				int a, b;
				if(i&3) {
					a = myRandInt(now.size());
					b = myRandInt(now.size()-1);
					if(a<=b) {
						++b;
					}
				}
				else if(i&4) {
					a = pena3.bottom();
					b = pena3.top();
				}
				else {
					a = pena2.bottom();
					b = pena2.top();
				}
				double nowA = now[a];
				double nowB = now[b];
				uint c = myRand();
				if(nowA<minNow[a] || nowB<minNow[b] || (c&1)) {
					if(minNow[a]<nowA && (c&2)) {
						double d = myRand2Float() * (nowA-minNow[a]) * paramA;
						nowA -= d;
						nowB += d;
					}
					else {
						if(nowB<minNow[b]) {
							continue;
						}
						double d = myRand2Float() * (nowB-minNow[b]) * paramA;
						nowB -= d;
						nowA += d;
					}
				}
				double aC = C[a] * pow(nowA, -1.0/eta);
				double bC = C[b] * pow(nowB, -1.0/eta);
				double aK = W[a] - aC;
				double bK = W[b] - bC;
				double aC2base = aC * pow(beta*(1-delta+alpha*Z2[a]*pow(aK, alpha-1.0)), 1.0/eta);
				double bC2base = bC * pow(beta*(1-delta+alpha*Z2[b]*pow(bK, alpha-1.0)), 1.0/eta);
				double aW2 = (1-delta) * aK + Z2[a] * pow(aK, alpha);
				double bW2 = (1-delta) * bK + Z2[b] * pow(bK, alpha);
				double diffBaseA = (log(aW2 * targetRate) - log(aC2base))*eta;
				double diffBaseB = (log(bW2 * targetRate) - log(bC2base))*eta;
				if(delta<0.2) {
					diffBaseA = -diffBaseA;
					diffBaseB = -diffBaseB;
				}
				double nextA = next[a];
				double nextB = next[b];
				double nextAB = nextA + nextB;
				double nextMinA = pow(aC2base / aW2, eta);
				double nextMinB = pow(bC2base / bW2, eta);

				g_diffBaseA = diffBaseA;
				g_diffBaseB = diffBaseB;
				g_nextAB = nextAB;
				nextA = searchConvex<false>(funcFast, nextMinA, nextAB-nextMinB);

				nextB = nextAB - nextA;

				double diffA = diffBaseA - log(nextA);
				double diffB = diffBaseB - log(nextB);
				double penaA = diffA*diffA;
				double penaB = diffB*diffB;

				if(penaA+penaB<pena[a]+pena[b]) {
					now[a] = nowA;
					now[b] = nowB;
					next[a] = nextA;
					next[b] = nextB;
					pena[a] = penaA;
					pena[b] = penaB;
					pena2.set(a, diffA);
					pena2.set(b, diffB);
					pena3.set(a, minNow[a]-nowA);
					pena3.set(b, minNow[b]-nowB);
				}
			}
			for(int i=0; i<(int)C.size(); ++i) {
				C[i] *= pow(now[i], -1.0/eta);
			}

			for(int i=0; i<(int)C.size(); ++i) {
				if(W[i]<=C[i]) {
					C[i] = (W[i] / (1+T-turn));
				}
			}
		}
		++turn;
		B = C;
		return C;
	}
};
