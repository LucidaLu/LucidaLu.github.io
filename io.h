#ifndef LUCIDA_IO
#define LUCIDA_IO
#include "universal.h"
NAMESPACE_BEGIN(io)
#ifdef OI
const int IN=1e6;
char in[IN],*ip=in,*ie=in;
#define getchar() ( ip==ie && (ie=(ip=in)+fread(in,1,IN,stdin),ip==ie)?EOF:*ip++ )
#endif
struct Istream {
	template <typename T>
	Istream &operator >>(T &);
	template <typename T>
	Istream &operator >>(T *);
}is;
template <typename T>
Istream &Istream::operator >>(T &x) {
	static bool neg;
	static char ch;
	neg=0;
	while(ch=getchar(),ch<'0' || '9'<ch) {
		neg|=ch=='-';
	}
	for(x=0;'0'<=ch && ch<='9';ch=getchar()) {
		(x*=10)+=ch-'0';
	}
	x=neg?-x:x;
	return *this;
}
template <>
Istream &Istream::operator >>(double &x) {
#ifndef OI
	scanf("%lf",&x);
#else
	throw std::invalid_argument("integer or char required");
#endif
	return *this;
}
template <>
Istream &Istream::operator >>(char &ch) {
	while(ch=getchar(),(ch=='\r' || ch=='\n' || ch==' ' || ch==EOF));
	return *this;
}
template <typename T>
Istream &Istream::operator >>(T *x) {
	char ch;
	while(ch=getchar(),(ch=='\r' || ch=='\n' || ch==' '  || ch==EOF));
	do *x++=ch;
	while(ch=getchar(),(ch!='\r' && ch!='\n' && ch!=' ' && ch!=EOF));
	*x=0;
	return *this;
}

#ifdef OO
const int OUT=1e6;
char out[OUT],*op=out,*oe=out+OUT;
#define flush() fwrite(out,1,op-out,stdout)
#define putchar(x) ( (op==oe?(flush(),op=out,*op++):*op++)=(x) )
#endif
struct Ostream {
#ifndef OO
	char dFormat[20];
	Ostream() {
		sprintf(dFormat,"%%lf");
	}
	void precision(const int &x) {
		sprintf(dFormat,"%%.%dlf",x);
	}
#else
	~Ostream() {
		flush();
	}
#endif
	template <typename T>
	Ostream &operator <<(T x);
	template <typename T>
	Ostream &operator <<(T *x);
}os;
template <typename T>
Ostream &Ostream::operator <<(T x) {
	if(x==0) {
		putchar('0');
	} else {
		if(x<0) {
			putchar('-');
			x=-x;
		}
		static char stack[100];
		static int top;
		for(top=0;x;x/=10) {
			stack[++top]=x%10+'0';
		}
		while(top) {
			putchar(stack[top--]);
		}
	}
	return *this;
}
template <>
Ostream &Ostream::operator <<(double x) {
#ifndef OO
	printf(dFormat,x);
#else
	throw std::invalid_argument("integer or char required");
#endif
	return *this;
}
template <>
Ostream &Ostream::operator <<(char x) {
	putchar(x);
	return *this;
}
template <typename T>
Ostream &Ostream::operator <<(T *x) {
	while(*x) {
		putchar(*x++);
	}
	return *this;
}
NAMESPACE_END
using io::is;
using io::os;
#endif
