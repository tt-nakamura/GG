// uses NTL
//   http://www.shoup.net/ntl

#ifndef __GG_h__
#define __GG_h__

#include<NTL/ZZ.h>

struct GG {// Gaussian integer x+iy
    NTL::ZZ x,y;// real and imaginary part
    GG() {;}
    GG(const NTL::ZZ& a, const NTL::ZZ& b) : x(a),y(b) {;}// a+bi
    GG(long a, long b) : x(a),y(b) {;}// a+bi
    GG(const NTL::ZZ& a) : x(a) {;}// a+0i
    GG(long a) : x(a) {;}// a+0i
    GG& operator=(const NTL::ZZ& a); // a+0i
    GG& operator=(long a);// a+0i
};

inline NTL::ZZ& real(GG& a) { return a.x; }
inline NTL::ZZ& imag(GG& a) { return a.y; }
inline const NTL::ZZ& real(const GG& a) { return a.x; }
inline const NTL::ZZ& imag(const GG& a) { return a.y; }

std::ostream& operator<<(std::ostream& s, const GG& a);
// print a as [a.x a.y]

void conv(GG& b, const NTL::ZZ& a);// b=a+0i
void conv(GG& b, long a);// b=a+0i

GG& operator+=(GG& b, const GG& a);// b+=a
GG& operator-=(GG& b, const GG& a);// b-=a
GG& operator*=(GG& b, const GG& a);// b*=a
GG& operator/=(GG& b, const GG& a);// b/=a
GG& operator%=(GG& b, const GG& a);// b%=a

inline long operator==(const GG& a, const GG& b) { return a.x==b.x && a.y==b.y; }
inline long operator!=(const GG& a, const GG& b) { return a.x!=b.x || a.y!=b.y; }
inline long operator==(const GG& a, long b) { return a.x==b && IsZero(a.y); }
inline long operator!=(const GG& a, long b) { return a.x!=b || !IsZero(a.y); }
inline long IsZero(const GG& a) { return NTL::IsZero(a.x) && NTL::IsZero(a.y); }
inline long IsOne(const GG& a) { return NTL::IsOne(a.x) && NTL::IsZero(a.y); }
inline long IsReal(const GG& a) { return NTL::IsZero(a.y); }

void set(GG& a); // a=1
void set(GG& a, const NTL::ZZ& x, const NTL::ZZ& y);// a=x+iy
void set(GG& a, long x, long y);// a=x+iy
void clear(GG& a);// a=0
void conj(GG& b, const GG& a);// b = x-iy when a==x+iy
void mirror(GG& b, const GG& a);// b = y+ix
void norm(NTL::ZZ& n, const GG& a);// n = x**2 + y**2
void negate(GG& b, const GG& a);// b = -a
void mul_i(GG& b, const GG& a);// b = a*i
void mul_i(GG& b, const GG& a, long k);// b = a*i^k
void div_i(GG& b, const GG& a);// b = a/i

long quadrant(const GG& a);
// return -1 if a==0
//         0 if real(a)>0 and imag(a)>=0
//         1 if real(a)<=0 and imag(a)>0
//         2 if real(a)<0 and imag(a)<=0
//         3 if real(a)>=0 and imag(a)<0

long FirstQuad(GG& b, const GG& a);
// b = a * i^k (k=0,1,2,3) such that Re(b)>0 and Im(b)>=0
// return k; k==0 if a==0

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

long divide2(GG& q, const GG& a);
// if a is divisible by 1+i, set q=a/(1+i) and return 1
// else return 0 (and q is unchanged)

void power(GG& b, const GG& a, long n);
// b = a^n; assume n>=0

long IsUnit(const GG& a);// test if |a|==1

long IsAssoc(const GG& a, const GG& b);
// return 1 if a = b*i^k for some k (k=0,1,2,3)
// else return 0

void GCD(GG&, const GG&, const GG&);
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
// by euclidean algorithm

void XGCD(GG&, GG&, GG&, const GG&, const GG&);
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
//   and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm

// a = random gaussian integer in square [-L,L]^2
void RandomBits(GG& a, long l);// 0 <= L < 2^l
void RandomLen(GG& a, long l);// 2^{l-1} <= L < 2^l
void RandomBnd(GG& a, const NTL::ZZ& n);// 0 <= L < n

void PowerMod(GG& b, const GG& a, long n, const GG& m);
void PowerMod(GG& b, const GG& a, const NTL::ZZ& n, const GG& m);
// b = a^n mod m; assume n>=0 and |a| < |m|

long ProbPrime(const GG& a, long NTRY=10);
// test if a is gaussian prime, i.e., return 1 if either
//   |a| is prime and |a|==3 (mod 4) or
//   |a|^2 is prime and |a|^2==1 (mod 4) or
//   |a|^2 == 2
// else return 0
// NTRY: number of trials of Miller-Rabin test

void GenPrime(GG& p, long l, long f=1, long err=80);
// generate random gaussian prime p.
// f must be 1 or 2; l must be l>=2
// if f==1, |p|^2=q and q==1 (mod 4)
// if f==2, p=q+0i and q==3 (mod 4)
// where q is random prime and 2^{l-1} <= q < 2^l
// probability of error is less than 2^-err
// p=x+iy is primary (x odd, y even, x+y==1 (mod 4))

long primary(GG& b, const GG& a);
// b = unit * a = x+iy such that
//   x==1 and y==0 or x==3 and y==2 (mod 4)
// return k such that a = i^k * b (k=0,1,2,3)
// Assume |a|^2 is odd

void ResSymb(GG& s, const GG& a, const GG& b);
// s = biquadratic residue symbol (a/b)_4 = 0,1,i,-1,-i
// Assume |a| < |b| and |b|^2 is odd
// Assume b is primary, but may not be prime

void FactorPrime(GG& a, const NTL::ZZ& p);
// given a prime number p where p==1 (mod 4),
// find x,y such that x^2 + y^2 = p
//   by cornacchia algorithm
// and return a = primary(x+yi)

void QrtRootMod(GG& x, const GG& a, const GG& p);
// solve x^4 == a (mod p)
// Assume p is primary prime and (a/p)_4 == 1
// Assume either
//   norm(p) is prime and norm(p)==1 (mod 4)
//   or norm(p) is prime^2 and p==3 (mod 4)

#endif // __GG_h__
