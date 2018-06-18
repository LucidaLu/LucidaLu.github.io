#ifndef LUCIDA_UNIVERSAL
#define LUCIDA_UNIVERSAL

#include <bits/stdc++.h>
/*
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <climits>
*/

#define array std::vector
#define NAMESPACE_BEGIN(X) namespace X {
#define NAMESPACE_END }

typedef long long int64;
typedef unsigned long long qword;
typedef long double decimal;

template <typename T>
inline bool chkmx(T &a,const T &b) {
	return a<b?(a=b,1):0;
}
template <typename T>
inline bool chkmn(T &a,const T &b) {
	return b<a?(a=b,1):0;
}
NAMESPACE_BEGIN(operator_ext)
template <typename T,typename C>
inline T &operator +=(T &a,const C &b) {
	return a=a+b;
}
template <typename T,typename C>
inline T &operator -=(T &a,const C &b) {
	return a=a-b;
}
template <typename T,typename C>
inline T &operator *=(T &a,const C &b) {
	return a=a*b;
}
template <typename T,typename C>
inline T &operator /=(T &a,const C &b) {
	return a=a/b;
}
NAMESPACE_END
const double eps=1e-6,PI=acos(-1.0);

inline int fcmp(const double &x) {
	return -eps<x && x<eps?0:(x<0?-1:1);
}

//initializer_list 比构造函数快
//参数 const int & int的空间应该是相同的，一样快。double 就不行了 
#endif
