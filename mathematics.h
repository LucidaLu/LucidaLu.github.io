#ifndef LUCIDA_MATHEMATICS
#define LUCIDA_MATHEMATICS
#include "universal.h"
#include "data_structure.h"
NAMESPACE_BEGIN(number_theory)
void Sift(int N,int prime[],int pcnt,int phi[],int mu[]) {//HDU 2874
	bool *notp=(bool*)calloc(N+1,sizeof(bool));
	for(int i=2;i<=N;++i) {
		if(!notp[i]) {
			prime[++pcnt]=i;
			phi[i]=i-1;
			mu[i]=-1;
		}
		for(int j=1;j<=pcnt && (int64)i*prime[j]<=N;++j) {
			int num=i*prime[j];
			notp[num]=1;
			if(i%prime[j]==0) {
				phi[num]=phi[i]*prime[j];
				mu[num]=0;
			} else {
				phi[num]=phi[i]*(prime[j]-1);
				mu[num]=-mu[i];
			}
		}
	}
	free(notp);
}
int Phi(int x) {//HDU 1286
	int mid=sqrt(x+0.5)+1,res=x;
	for(int i=2;i<=mid;++i) if(x%i==0) {
		(res/=i)*=(i-1);
		do x/=i;
		while(x%i==0);
	}
	if(x!=1) {
		(res/=x)*=(x-1);
	}
	return res;
}
int gcd(int a,int b) {
	for(int t;b;t=a,a=b,b=t%b);
	return a;
}
int Exgcd(int a,int b,int64 &x,int64 &y) {
	if(!b) {
		x=1,y=0;
		return a;
	} else {
		int d=Exgcd(b,a%b,x,y);int64 t;
		t=x;x=y;y=t-(a/b)*y;
		return d;
	}
}
int ModEquation(int a,int b,int P) {//ax\equiv b\mod P BZOJ 2242
	int64 x,y;
	int d=Exgcd(a,P,x,y);
	if(b%d) {
		return -1;
	} else {
		x*=b/d;int k=P/d;
		((x%=k)+=k)%=k;
		return x;
	}
}
int Pow(int64 base,int index,int P) {
	int64 res;
	for(res=1;index;(base*=base)%=P,index>>=1) {
		if(index&1) {
			(res*=base)%=P;
		}
	}
	return res;
}
int Inverse(int a,int P) {
	int64 x,y;
	int d=Exgcd(a,P,x,y);
	//assert(d==1);
	if(d!=1) {
		throw;
	}
	return (x%P+P)%P;
}
int Lucas(int n,int m,int P) {//HDU 3037
	if(n<m) {
		return 0;
	}
	int *fact=(int*)malloc((P+1)*sizeof(int));
	fact[0]=1;
	for(int i=1;i<=P;++i) {
		fact[i]=(int64)fact[i-1]*i%P;
	}
	int64 res=1;
	for(;n && m;n/=P,m/=P) {
		int a=n%P,b=m%P;
		if(a<b) {
			res=0;
			break;
		}
		(res*=(int64)fact[a]*Inverse(fact[b],P)%P*Inverse(fact[a-b],P)%P)%=P;
	}
	free(fact);
	return res;
}
int CRT(const array<int> &A,const array<int> &m,int M) {
	//A[i]*t[i]*(M/m[i])
	int64 res=0;
	for(int i=0,sz=A.size();i<sz;++i) {
		int Mi=M/m[i];
		(res+=(int64) A[i] * Inverse(Mi,m[i])%M * Mi %M )%=M;
	}
	return res;
}

NAMESPACE_BEGIN(exlucas)
int *partFact,p,t,pt;
int Calc(int n) {
	return n<p?partFact[n]:(int64)Pow(partFact[pt],n/pt,pt)*partFact[n%pt]%pt*Calc(n/p)%pt;	
}
int C(const int n,const int m) {
	int index=0;
	for(int x=n;x;index+=(x=x/p));
	for(int x=m;x;index-=(x=x/p));
	for(int x=n-m;x;index-=(x=x/p));
	if(index>=t) {
		return 0;
	} else {
		partFact[0]=1;
		for(int i=1;i<=pt;++i) {
			partFact[i]=i%p==0?partFact[i-1]:(int64)partFact[i-1]*i%pt;
		}
		int fn=Calc(n),fm=Calc(m),fnm=Calc(n-m);
		return (int64)Pow(p,index,pt)*fn%pt*Inverse(fm,pt)%pt*Inverse(fnm,pt)%pt;
	}
}
int ExLucas(int n,int m,const int P) {//BZOJ 3283
	if(n<m) {
		return 0;
	}
	array<int> A,dpt,dp,dt;
	int x=P,mid=sqrt(x+0.5)+1,maxpt=0;
	for(int i=2;i<=x;i=i==mid?(x==1?INT_MAX:x):i+1) if(x%i==0) {
		p=i,t=0,pt=1;
		do x/=i,t++,pt*=i;
		while(x%i==0);
		dp.push_back(p);
		dt.push_back(t);
		dpt.push_back(pt);
		chkmx(maxpt,pt);
	}
	partFact=(int*)malloc((maxpt+1)*sizeof(int));
	A.resize(dp.size());
	for(int i=0,sz=dp.size();i<sz;++i) {
		p=dp[i];t=dt[i];pt=dpt[i];
		A[i]=C(n,m);
	}
	free(partFact);
	return CRT(A,dpt,P);
}
NAMESPACE_END
using exlucas::ExLucas;
int BSGS(int a,int b,int P) {//a^x\equiv b \mod P BZOJ 2242
	a%=P;b%=P;
	if(a==0) {
		return -1;
	}
	data_structure::HashMap S;
	int64 baj=b;
	int mid=sqrt(P+0.5)+1;
	for(int j=0;j<=mid;++j,(baj*=a)%=P) {
		S.Insert(baj,j);
	}
	int64 asp=Pow(a,mid,P),aspi=asp;
	for(int i=1;i<=mid;++i,(aspi*=asp)%=P) {
		int p=S.Find(aspi);
		if(p!=-1) {
			return i*mid-p;
		}
	}
	return -1;
}
int ExBSGS(int a,int b,int P) {//BZOJ 3283
	a%=P;b%=P;
	if(a==0) {//一定要加特判！ 0^0 is undefined! 
		return -1;
	} else if(b==1) {
		return 0;
	}
	int t=0;int64 k=1;
	for(int d;(d=gcd(a,P))!=1;) {
		if(b%d) {
			return -1;
		} else {
			b/=d;P/=d;
			(k*=a/d)%=P;
			++t;
			if(k==b) {
				return t;
			}
		}
	}
	data_structure::HashMap S;
	int64 baj=b;
	int mid=sqrt(P+0.5)+1;
	for(int j=0;j<=mid;++j,(baj*=a)%=P) {
		S.Insert(baj,j);//md 存反了 
	}
	int64 asp=Pow(a,mid,P),kaspi=k*asp%P;
	for(int i=1;i<=mid;++i,(kaspi*=asp)%=P) {
		int j=S.Find(kaspi);
		if(j!=-1) {
			return mid*i-j+t;
		}
	}
	return -1;
}
void Divide(int x,int *p,int *t,int &cnt) {
	int mid=sqrt(x+0.5)+1;cnt=0;
	for(int i=2;i<=mid;++i) if(x%i==0) {
		p[++cnt]=i;t[cnt]=0;
		do x/=i,t[cnt]++;
		while(x%i==0);
	}
	if(x!=1) {
		p[++cnt]=x;t[cnt]=1;
	}
}
int PrimitiveRoot(int P) {//51Nod 1135
	int sqrtP=sqrt(P+0.5)+1;
	int *p=(int*)malloc((sqrtP)*sizeof(int)),*t=(int*)malloc((sqrtP)*sizeof(int)),cnt,Ans=-1,phiP=Phi(P);
	Divide(phiP,p,t,cnt);
	for(int g=2;g<P;++g) {
		bool tak=1;
		for(int i=1;i<=cnt;++i) {
			if(Pow(g,phiP/p[i],P)==1) {
				tak=0;
				break;
			}
		}
		if(tak) {
			Ans=g;
			break;
		}
	}
	free(p);free(t);
	return Ans;
}
NAMESPACE_END

NAMESPACE_BEGIN(linear_algebra)
bool isPositive(double const &x) {
	return fcmp(x)>0;
}
bool isNegative(double const &x) {
	return fcmp(x)<0;
}
struct Simplex {//在uoj上被卡精度了..eps=1e-7才可以. 
	static const double IMPOSSIBLE,UNLIMITED;
	double **a,*b,*c,Ans;
	int n,m,cnt,*idvar,*idine;
	Simplex(int n,int m):n(n),m(m),cnt(0) {
		Ans=0;
		c=(double*)calloc(n+1,sizeof(double));
		a=(double**)calloc(m+1,sizeof(double*));
		for(int i=1;i<=m;++i) {
			a[i]=(double*)calloc(n+1,sizeof(double));
		}
		b=(double*)calloc(m+1,sizeof(double));
		idvar=(int*)calloc(n+1,sizeof(int));
		idine=(int*)calloc(m+1,sizeof(int));
	}
	~Simplex() {
		free(c);
		for(int i=1;i<=m;++i) {
			free(a[i]);
		}
		free(a);free(b);
		free(idvar);free(idine);
	}
	void Set(double nc[]) {
		std::copy(nc+1,nc+1+n,c+1);
		for(int i=1;i<=n;++i) {
			idvar[i]=i;
		}
	}
	void Add(double na[],double nb) {
		if(++cnt>m) {
			throw;
		}
		std::copy(na+1,na+1+n,a[cnt]+1);
		b[cnt]=nb;
		idine[cnt]=n+cnt;
	}
	bool Init() {
		int l,e;
		srand(0x1f1f1f1f);//一定要加随机扰动 
		while((l=std::find_if(b+1,b+1+m,isNegative)-b)!=m+1) {
			if((e=std::find_if(a[l]+1,a[l]+1+n,isNegative)-a[l])!=n+1) {
				for(int i=1;i<=m;++i) {
					if(isNegative(b[i]) && rand()&1) {
						l=i;
					}
				}
				for(int j=1;j<=n;++j) {
					if(isNegative(a[l][j]) && rand()&1) {
						e=j;
					}
				}
				Pivot(l,e);
			} else {
				return 0;
			}
		}
		return 1;
	}
	void Pivot(int l,int e) {//l是不等式序号 e是变量序号 
		std::swap(idine[l],idvar[e]);
		for(int j=1;j<=n;++j) if(j!=e) {
			a[l][j]/=a[l][e];
		}
		b[l]/=a[l][e];
		a[l][e]=1/a[l][e];
		for(int i=1;i<=m;++i) if(i!=l && fcmp(a[i][e])) {
			b[i]-=a[i][e]*b[l];
			for(int j=1;j<=n;++j) if(j!=e) {
				a[i][j]-=a[i][e]*a[l][j];
			}
			a[i][e]=-a[i][e]*a[l][e];
		}
		Ans+=c[e]*b[l];
		for(int j=1;j<=n;++j) if(j!=e) {
			c[j]-=a[l][j]*c[e];
		}
		c[e]=-c[e]*a[l][e];//*不是/ 
	}
	double Optimum(double x[]) {
		if(!Init()) {
			return IMPOSSIBLE;
		}
		int l,e;
		while((e=std::find_if(c+1,c+1+n,isPositive)-c)!=n+1) {
			double lim=UNLIMITED;
			for(int i=1;i<=m;++i) {
				if(isPositive(a[i][e]) && chkmn(lim,b[i]/a[i][e])) {
					l=i;
				}
			}
			if(lim==UNLIMITED) {
				return UNLIMITED;
			} else {
				Pivot(l,e);
			}
		}
		std::fill(x+1,x+1+n,0);
		for(int i=1;i<=m;++i) {
			if(idine[i]<=n) {
				x[idine[i]]=b[i];
			}
		}
		return Ans;
	}
};
const double Simplex::IMPOSSIBLE=DBL_MIN,Simplex::UNLIMITED=DBL_MAX; 
NAMESPACE_END

NAMESPACE_BEGIN(polynomial)
template <typename T>
void Rader(T a[],int n) {
	for(int i=1,j=n>>1;i<n-1;++i) {
		if(i<j) {
			std::swap(a[i],a[j]);
		}
		int k=n>>1;
		while(j>=k) {
			j-=k;
			k>>=1;
		}
		j+=k;
	}
}
#ifdef DFT_IMPL
struct complex {
	double r,i;
	friend complex operator +(const complex &a,const complex &b) {
		return (complex){a.r+b.r,a.i+b.i};
	}
	friend complex operator -(const complex &a,const complex &b) {
		return (complex){a.r-b.r,a.i-b.i};
	}
	friend complex operator *(const complex &a,const complex &b) {
		return (complex){a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};
	}
	friend complex operator /(const complex &a,const double &b) {
		return (complex){a.r/b,a.i/b};
	}
};
void DFT(complex a[],int n,int op) {
	using namespace operator_ext;
	Rader(a,n);
	for(int len=2;len<=n;len<<=1) {
		complex rt={cos(2*PI*op/len),sin(2*PI*op/len)};
		for(int i=0;i<n;i+=len) {
			complex w={1,0};
			for(int j=i;j<i+(len>>1);++j) {
				complex u=a[j],t=w*a[j+(len>>1)];
				a[j]=u+t,a[j+(len>>1)]=u-t;
				w*=rt;
			}
		}
	}
	if(op==-1) {
		for(int i=0;i<n;++i) {
			a[i]/=n;
		}
	}
}
#else
struct NumberTheoricTransform {
	int P,g,index2;int64 *wn;
	NumberTheoricTransform(int P):P(P) {
		g=number_theory::PrimitiveRoot(P);
		index2=0;
		for(int x=P-1;~x&1;x>>=1) {
			index2++;
		}
		wn=(int64*)malloc((index2+1)*sizeof(int64));
		wn[0]=1;
		for(int i=1;i<=index2;++i) {
			wn[i]=number_theory::Pow(g,(P-1)/(1<<i),P);
		}
	}
	~NumberTheoricTransform() {
		free(wn);
	}
	void operator () (int64 *a,int n,int op) {
		if((1<<index2)<n) {
			throw;
		}
		Rader(a,n);
		for(int len=2,index=1;len<=n;len<<=1,++index) {//index..
			int64 rt=wn[index];//不要写成了len.. 
			for(int i=0;i<n;i+=len) {
				int64 w=1;
				for(int j=i;j<i+(len>>1);++j) {
					int64 u=a[j],t=w*a[j+(len>>1)]%P;//不要忘记*w... 
					a[j]=(u+t)%P;a[j+(len>>1)]=(u-t)%P;//前+后-不要写反了。。 
					(w*=rt)%=P;
				}
			}
		}
		if(op==-1) {
			std::reverse(a+1,a+n);
			int invN=number_theory::Inverse(n,P);
			for(int i=0;i<n;++i) {
				(a[i]*=invN)%=P;
			}
		}
	}
};
#endif
struct Polynomial {
	int *a,n;//n是次数界 
	Polynomial(int n=0):a((int*)calloc(n,sizeof(int))),n(n) {}
	Polynomial(const Polynomial &A) {
		*this=A;
	}
	~Polynomial() {
		free(a);
	}
	Polynomial &operator =(const Polynomial &A) {
		n=A.n;
		a=(int*)malloc(n*sizeof(int));
		memcpy(a,A.a,n*sizeof(int));
		return *this;
	}
	int &operator [](const int &x) {
		return a[x];
	}
	const int &operator [](const int &x) const {
		return a[x];
	}
	friend Polynomial operator *(const Polynomial &A,const Polynomial &B) {
		using namespace operator_ext;
		int n=1;
		Polynomial C(A.n+B.n-1);
		while(n<C.n) {
			n<<=1;
		}
#ifdef DFT_IMPL
		complex *a=(complex*)calloc(n,sizeof(complex)),*b=(complex*)calloc(n,sizeof(complex));
		for(int i=0;i<n;++i) {

			a[i]=(complex){(i<A.n?A[i]:0.0),0.0};
			b[i]=(complex){(i<B.n?B[i]:0.0),0.0};
		}
		DFT(a,n,1);DFT(b,n,1);
		for(int i=0;i<n;++i) {
			a[i]*=b[i];
		}
		DFT(a,n,-1);
		for(int i=0;i<C.n;++i) {
			C[i]=(int)(a[i].r+0.5);
		}
#else
		static const int P=104857601; 
		NumberTheoricTransform NTT(P);
		int64 *a=(int64*)calloc(n,sizeof(int64)),*b=(int64*)calloc(n,sizeof(int64));
		for(int i=0;i<n;++i) {
			a[i]=i<A.n?A[i]:0;
			b[i]=i<B.n?B[i]:0;
		}
		NTT(a,n,1);NTT(b,n,1);
		for(int i=0;i<n;++i) {
			(a[i]*=b[i])%=P;
		}
		NTT(a,n,-1);
		for(int i=0;i<C.n;++i) {
			C[i]=(a[i]+P)%P;
		}
#endif
		free(a);free(b);
		return C;
	}
};

NAMESPACE_END
NAMESPACE_BEGIN(combinatorics)
using number_theory::ExLucas;
//Polya
//Burnside
NAMESPACE_END

NAMESPACE_BEGIN(game_theory)
//SG
NAMESPACE_END

NAMESPACE_BEGIN(graph_theory)
struct Edge {
	int to,v;Edge *pre;
	Edge(int to,int v,Edge *pre):to(to),v(v),pre(pre) {}
};
struct Graph {
	typedef unsigned int val_t;
	int n;Edge **G;
	Graph(int n):n(n),G(new Edge*[n+1]()) {}
	~Graph() {
		for(int i=1;i<=n;++i) {
			for(Edge *pre,*e=G[i];e;e=pre) {
				pre=e->pre;
				delete e;
			}
		}
		delete[] G;
	}
	Edge *operator [](const int &x) {
		return G[x];
	}
};//还有弦图！！ 
struct DiGraph : Graph {
	DiGraph(int n):Graph(n) {}
	void operator ()(int f,int t,val_t v=1) {
		G[f]=new Edge(t,v,G[f]);
	}
	void SCC(array<array<int> > &);
};
struct UndiGraph : Graph {
	UndiGraph(int n):Graph(n) {}
	void operator ()(int f,int t,val_t v=1) {
		G[f]=new Edge(t,v,G[f]);
		G[t]=new Edge(f,v,G[t]);
	}
};
struct TarjanSCC {
	DiGraph &G; 
	int *&sn,&sc,*stack,top,*low,*dfn,dc;
	bool *instack;
	TarjanSCC(DiGraph &G,int *sn,int &sc):G(G),sn(sn),sc(sc) {
		top=dc=sc=0;
		stack=new int[G.n+1]();
		low=new int[G.n+1];
		dfn=new int[G.n+1]();
		instack=new bool[G.n+1](); 
		for(int i=1;i<=G.n;++i) if(!dfn[i]) {
			DFS(i,0);
		}
		delete[] stack;
		delete[] low;
		delete[] dfn;
		delete[] instack;
	}
	void DFS(int pos,int fa) {
		stack[++top]=pos;
		low[pos]=dfn[pos]=++dc;
		instack[pos]=1;
		for(Edge *e=G[pos];e;e=e->pre) if(e->to!=fa) {
			int u=e->to;
			if(!dfn[u]) {
				DFS(u,pos);
				chkmn(low[pos],low[u]);
			} else if(instack[u]) {
				chkmn(low[pos],dfn[u]);
			}
		}
		if(low[pos]==dfn[pos]) {
			sc++;
			while(stack[top+1]!=pos) {
				sn[stack[top]]=sc;
				instack[stack[top]]=0;
				top--;
			}
		}
	}
};
void DiGraph::SCC(array<array<int> > &scc) {
	int *sn=new int[n+1](),sc=0;
	TarjanSCC(*this,sn,sc);
	scc.resize(sc+1);
	for(int i=1;i<=n;++i) {
		scc[sn[i]].push_back(i);
	}
	delete[] sn;
}
NAMESPACE_END
#endif
