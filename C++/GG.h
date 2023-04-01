#ifndef __GG_h__
#define __GG_h__

#include<NTL/ZZ.h>
using namespace NTL;

struct GG {// Gaussian integer x+iy
    ZZ x,y;
    GG() {;}
    GG(const GG& a) : x(a.x), y(a.y) {;}
    GG(const ZZ& a, const ZZ& b) : x(a),y(b) {;}
    GG(long a, long b) : x(a),y(b) {;}
    GG(const ZZ& a) : x(a) {;}
    GG(long a) : x(a) {;}
    GG& operator=(const GG&);
    GG& operator=(const ZZ&);
    GG& operator=(long);
};

inline ZZ& real(GG& a) { return a.x; }
inline ZZ& imag(GG& a) { return a.y; }
inline const ZZ& real(const GG& a) { return a.x; }
inline const ZZ& imag(const GG& a) { return a.y; }

const GG& operator+=(GG&, const GG&);
const GG& operator-=(GG&, const GG&);
const GG& operator*=(GG&, const GG&);
const GG& operator/=(GG&, const GG&);
const GG& operator%=(GG&, const GG&);

inline long operator==(const GG& a, const GG& b) { return a.x==b.x && a.y==b.y; }
inline long IsZero(const GG& a) { return IsZero(a.x) && IsZero(a.y); }
inline long IsOne(const GG& a) { return IsOne(a.x) && IsZero(a.y); }
inline long IsReal(const GG& a) { return IsZero(a.y); }

void set(GG&);
void set(GG&, const ZZ&, const ZZ&);
void set(GG&, long, long);
void clear(GG& a);
void conj(GG&, const GG&);
void mirror(GG&, const GG&);
void norm(ZZ&, const GG&);
void negate(GG&, const GG&);
void mul_i(GG&, const GG&);
void mul_i(GG&, const GG&, long);
void div_i(GG&, const GG&);

long quadrant(const GG&);
long FirstQuad(GG&, const GG&);

void add(GG&, const GG&, const GG&);
void sub(GG&, const GG&, const GG&);
void mul(GG&, const GG&, const GG&);
void sqr(GG&, const GG&);
void div(GG&, const GG&, const GG&);
void rem(GG&, const GG&, const GG&);
void DivRem(GG&, GG&, const GG&, const GG&);
long divide(GG&, const GG&, const GG&);
long divide(const GG&, const GG&);
void power(GG&, const GG&, long);

long IsUnit(const GG&);
long IsAssoc(const GG&, const GG&);
long ProbPrime(const GG&);

void RandomBnd(GG&, const ZZ&);
void RandomBits(GG&, long);
void RandomLen(GG&, long);
void GenPrime(GG&, long, long=80);
void GenRealPrime(GG&, long, long=80);

void GCD(GG&, const GG&, const GG&);
void XGCD(GG&, GG&, GG&, const GG&, const GG&);

std::ostream& operator<<(std::ostream&, const GG&);

void cornacchia(ZZ&, ZZ&, const ZZ&);

#endif // __GG_h__
