#ifndef LUCIDA_GEOMETRY2D
#define LUCIDA_GEOMETRY2D
#include "universal.h"
NAMESPACE_BEGIN(geometry2D)
struct point {
	double x,y;
	point() {}
	point(double x,double y):x(x),y(y) {}
	double norm() const {
		return sqrt(x*x+y*y);
	}
	point normal() const {
		return point(-y,x)/norm();
	}
	friend point operator +(const point &a,const point &b) {
		return point(a.x+b.x,a.y+b.y);
	}
	friend point operator -(const point &a,const point &b) {
		return point(a.x-b.x,a.y-b.y);
	}
	friend point operator *(const point &a,const double &k) {
		return point(a.x*k,a.y*k);
	}
	friend point operator /(const point &a,const double &k) {
		return point(a.x/k,a.y/k);
	}
	friend double Inner(const point &a,const point &b) {
		return a.x*b.x+a.y*b.y;
	}
	friend double Outer(const point &a,const point &b) {
		return a.x*b.y-a.y*b.x;
	}
};
typedef point vector;
typedef point* polygon;
struct line {
	point p;
	vector v;
	double ang;
	line() {}
	line(point p,vector v):p(p),v(v) {
		ang=atan2(v.y,v.x);
	}
};
typedef line* plane;
struct circle {
	point o;
	double r;
	circle() {}
	circle(point o,double r):o(o),r(r) {}
};
double Dist(const point &a,const point &b) {
	return (a-b).norm();
}
double Dist(const point &a,const line &l) {
	return fabs(Outer(a-l.p,l.v))/l.v.norm();
}
point Cross(const line &ln1,const line &ln2) {
	double t=Outer(ln2.v,ln1.p-ln2.p)/Outer(ln1.v,ln2.v);
	return ln1.p+ln1.v*t;
}

bool cmpXY(const point &a,const point &b) {
	return a.x!=b.x?a.x<b.x:a.y<b.y; 
}
void ConvexHull(polygon p,int n,polygon hull,int &hc) {//POJ 1113 	
	/*
	stack需要2倍空间 
	k=top+1
	转的方向 
	*/ 
	using std::sort;
	using std::copy;
	point *stack=(point*)malloc((n<<1)*sizeof(point));
	int top=0;
	sort(p+1,p+1+n,cmpXY);
	for(int i=1;i<=n;++i) {
		while(top>=2 && fcmp(Outer(stack[top]-stack[top-1],p[i]-stack[top]))<=0) {//sqrt -。。 但是这样不科学 
			top--;
		}
		stack[++top]=p[i];
	}
	int k=top+1;
	for(int i=n-1;i;--i) {
		while(top>=k && fcmp(Outer(stack[top]-stack[top-1],p[i]-stack[top]))<=0) {
			top--;
		}
		stack[++top]=p[i];
	}
	if(top!=1) {
		top--;
	}
	copy(stack+1,stack+1+top,hull+1);
	hc=top;
	free(stack);
}

double MaxDist(polygon p,int n) {//POJ 2187
	//输入必须有序
	//好吧 宏一定要在上下文定义！ 
	if(n==2) {//一定要判线的情况！ 
		return Dist(p[1],p[2]);
	} else {
		double res=0;
		#define succ(i) ((i)+1>n?1:(i)+1)
		for(int i=1,cp=2;i<=n;++i) {
			line cl(p[i],p[succ(i)]-p[i]);
			while(Dist(p[cp],cl)<Dist(p[succ(cp)],cl)) {
				cp=succ(cp);
			}
			chkmx(res,std::max(Dist(p[cp],p[i]),Dist(p[cp],p[succ(i)])));
		}
		#undef succ
		return res;
	}
}

bool cmpAngle(const line &a,const line &b) {
	return a.ang<b.ang;
}
bool OnLeft(const point &p,const line &ln) {
	return fcmp(Outer(ln.v,p-ln.p))>0;
}
bool Paral(const line &ln1,const line &ln2) {
	return fcmp(Outer(ln1.v,ln2.v))==0;
}
void Intersect(plane ln,int n,plane uni,int &uc) {//POJ 3525
	std::sort(ln+1,ln+1+n,cmpAngle);
	point *Qp=(point*)malloc((n+1)*sizeof(point));
	line *Qln=(line*)malloc((n+1)*sizeof(line));
	int he,ta;Qln[he=ta=1]=ln[1];
	#define count(a,b) ((b)-(a)+1)
	for(int i=2;i<=n;++i) {
		while(count(he,ta)>=2 && !OnLeft(Qp[ta-1],ln[i])) {
			ta--;
		}
		while(count(he,ta)>=2 && !OnLeft(Qp[he],ln[i])) {
			he++;
		}
		Qln[++ta]=ln[i];
		if(Paral(Qln[ta-1],Qln[ta])){
			--ta;
			if(OnLeft(ln[i].p,Qln[ta])) {
				Qln[ta]=ln[i];
			}
		}
		if(count(he,ta)>=2) {
			Qp[ta-1]=Cross(Qln[ta-1],Qln[ta]);
		}
	}
	while(count(he,ta)>=2 && !OnLeft(Qp[ta-1],Qln[he])) {
		ta--;
	}
	if(count(he,ta)>=2) {
		std::copy(Qln+he,Qln+ta+1,uni+1);
		uc=count(he,ta);
	} else {
		uc=0;
	}
}
NAMESPACE_END
#endif
