// uses NTL
//   http://www.shoup.net/ntl

#ifndef __GG_h__
#define __GG_h__

#include<NTL/ZZ.h>

struct GG {// Gaussian integer x+iy
    NTL::ZZ x,y;// real and imaginary part
    GG() {;}
    GG(const NTL::ZZ& a, const NTL::ZZ& b) : x(a),y(b) {;}
    GG(long a, long b) : x(a),y(b) {;}
    GG(const NTL::ZZ& a) : x(a) {;}// a+0i
    GG(long a) : x(a) {;}// a+0i
    GG& operator=(const NTL::ZZ& a); // a+0i
    GG& operator=(long a);// a+0i
};

inline NTL::ZZ& real(GG& a) { return a.x; }
inline NTL::ZZ& imag(GG& a) { return a.y; }
inline const NTL::ZZ& real(const GG& a) { return a.x; }
inline const NTL::ZZ& imag(const GG& a) { return a.y; }

GG& operator+=(GG& b, const GG& a);
GG& operator-=(GG& b, const GG& a);
GG& operator*=(GG& b, const GG& a);
GG& operator/=(GG& b, const GG& a);
GG& operator%=(GG& b, const GG& a);

inline long operator==(const GG& a, const GG& b) { return a.x==b.x && a.y==b.y; }
inline long IsZero(const GG& a) { return NTL::IsZero(a.x) && NTL::IsZero(a.y); }
inline long IsOne(const GG& a) { return NTL::IsOne(a.x) && NTL::IsZero(a.y); }
inline long IsReal(const GG& a) { return NTL::IsZero(a.y); }

void set(GG& a); // a=1
void set(GG& a, const NTL::ZZ& x, const NTL::ZZ& y);// a=x+iy
void set(GG& a, long x, long y);// a=x+iy
void clear(GG& a);// a=9
void conj(GG& b, const GG& a);// b = x-iy when a==x+iy
void mirror(GG& b, const GG& a);// b = y+ix
void norm(NTL::ZZ& n, const GG& a);// n = x**2 + y**2
void negate(GG& b, const GG& a);// b = -a
void mul_i(GG& b, const GG& a);// b = a*i
void mul_i(GG& b, const GG& a, long e);// b = a*i^e
void div_i(GG& b, const GG& a);// b = a/i

long quadrant(const GG& a);
// return -1 if a==0
//         0 if real(a)>0 and imag(a)>=0
//         1 if real(a)<=0 and imag(a)>0
//         2 if real(a)<0 and imag(a)<=0
//         3 if real(a)>=0 and imag(a)<0

long FirstQuad(GG& b, const GG& a);
// b = a * i^e (e=0,1,2,3) such that Re(b)>0 and Im(b)>=0
// return e; e==0 if a==0

void primary(GG& b, const GG& a);
// b = a*unit == x+yi to make
//   x odd and y even and x+y==1 (mod 4)
// assume norm(a) != 0 (mod 2)

void add(GG& c, const GG& a, const GG& b);// c=a+b
void sub(GG& c, const GG& a, const GG& b);// c=a-b
void mul(GG& c, const GG& a, const GG& b);// c=a*b
void sqr(GG& b, const GG& a);// b=a*a

void div(GG& q, const GG& a, const GG& b);
// q = quotient of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2

void rem(GG& r, const GG& a, const GG& b);
// r = remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2

void DivRem(GG& q, GG& r, const GG& a, const GG& b);
// q,r = quotient and remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2

long divide(GG& q, const GG& a, const GG& b);
// if a/b is divisible, set q=a/b and return 1
// else return 0 (and q is unchanged)

long divide(const GG& a, const GG& b);
// if a/b is divisible, return 1, else return 0

void power(GG&, const GG&, long);
// b = a^n; assume n>=0

long IsUnit(const GG& a);// test if a==1
long IsAssoc(const GG& a, const GG& b);// test if |a/b|==1

long ProbPrime(const GG& a, long NTRY=10);
// test if a is gaussian prime, i.e.,
//   |a|^2 is prime and |a|^2==1 (mod 4) or
//   Im(a)==0 and |Re(a)| is prime and |Re(a)|==3 (mod 4) or
//   Re(a)==0 and |Im(a)| is prime and |Im(a)|==3 (mod 4)
// NTRY: number of trials of Miller-Rabin test

void GCD(GG&, const GG&, const GG&);
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
// by euclidean algorithm

void XGCD(GG&, GG&, GG&, const GG&, const GG&);
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
//   and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm

void RandomBnd(GG&, const NTL::ZZ&);
// a = random gaussian integer such that
//   |Re(a)|<n and |Im(a)|<n

void GenPrime(GG&, long, long=1, long=80);
// generate random gaussian prime.
// p = gaussian integer such that |p|^2 == q^f
//   where q is random prime integer
// l = bit length of q, so that 2^{l-1} < q < 2^l
// f = 1 or 2 (P == 1 or 3 mod 4, respectively).
//   If f==1, p is imaginary, Re(p)>0, Im(p)>0.
//   If f==2, p = q (real and positive).
// err = bound for error probability < 2^{-err}
// If f<1 or f>2, f is set to 1 or 2, respectively.

std::ostream& operator<<(std::ostream& s, const GG& a);
// print a as [a.x a.y]

void FactorPrime(GG& a, const NTL::ZZ& p);
// given a prime number p where p==1 (mod 4),
// find x,y such that x^2 + y^2 = p
//   by cornacchia algorithm
// and return a = primary(x+yi)

#endif // __GG_h__
