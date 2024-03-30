// uses NTL
//   http://www.shoup.net/ntl

#include "GG.h"
using namespace NTL;

GG& GG::operator=(const ZZ& a)
{ x=a; clear(y); return *this; }// a+0i

GG& GG::operator=(long a)
{ x=a; clear(y); return *this; }// a+i0

GG& operator+=(GG& a, const GG& b) {// a+=b
    a.x += b.x;
    a.y += b.y;
    return a;
}

GG& operator-=(GG& a, const GG& b) {// a-=b
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

GG& operator*=(GG& a, const GG& b) {// a*=b
    mul(a,a,b);
    return a;
}

GG& operator/=(GG& a, const GG& b) {// a/=b
    div(a,a,b);
    return a;
}

GG& operator%=(GG& a, const GG& b) {// a%=b
    rem(a,a,b);
    return a;
}

void set(GG& a) { set(a.x); clear(a.y); }// a=1
void set(GG& a, const ZZ& x, const ZZ& y) { a.x=x; a.y=y; }// a=x+iy
void set(GG& a, long x, long y) { a.x=x; a.y=y; }// a=x+iy

void clear(GG& a) { clear(a.x); clear(a.y); }// a=0

void conj(GG& b, const GG& a) {// b = complex conjugate of a
    if(&b!=&a) b.x = a.x;
    negate(b.y, a.y);
}

void mirror(GG& b, const GG& a) {// b = Im(a) + i*Re(a)
    if(&b==&a) swap(b.x, b.y);
    else { b.x = a.y; b.y = a.x; }
}

void norm(ZZ& b, const GG& a) {// b = |a|^2
    ZZ s,t;
    sqr(s, a.x);
    sqr(t, a.y);
    add(b,s,t);
}

void negate(GG& b, const GG& a) {// b=-a
    negate(b.x, a.x);
    negate(b.y, a.y);
}

void mul_i(GG& b, const GG& a) {// b = i*a
    if(&b==&a) {
        swap(b.x, b.y);
        negate(b.x, b.x);
    }
    else {
        negate(b.x, a.y);
        b.y = a.x;
    }
}

void div_i(GG& b, const GG& a) {// b=a/i
    if(&b==&a) {
        swap(b.x, b.y);
        negate(b.y, b.y);
    }
    else {
        b.x = a.y;
        negate(b.y, a.x);
    }
}

void mul_i(GG& b, const GG& a, long e) {// b = a * i^e
    if((e&=3)==1) mul_i(b,a);
    else if(e==2) negate(b,a);
    else if(e==3) div_i(b,a);
    else if(&b!=&a) b=a;
}

long quadrant(const GG& a)
// return -1 if a==0
//         0 if real(a)>0 and imag(a)>=0
//         1 if real(a)<=0 and imag(a)>0
//         2 if real(a)<0 and imag(a)<=0
//         3 if real(a)>=0 and imag(a)<0
{
    if(IsZero(a.x)) {
        if(IsZero(a.y)) return -1;
        else if(sign(a.y) > 0) return 1;
        else return 3;
    }
    else if(sign(a.x) > 0) {
        if(sign(a.y) < 0) return 3;
        else return 0;
    }
    else if(sign(a.y) > 0) return 1;
    else return 2;
}

long FirstQuad(GG& b, const GG& a)
// b = a * i^e (e=0,1,2,3) such that Re(b)>0 and Im(b)>=0
// return e; e==0 if a==0
{
    long c(quadrant(a));
    if(c>0) c=4-c; else c=0;
    mul_i(b,a,c);
    return c;
}

void primary(GG& b, const GG& a)
// b = a*unit == x+yi to make
//   x odd and y even and x+y==1 (mod 4)
// assume norm(a) != 0 (mod 2)
{
    long x(a.x%4), y(a.y%4);
    if(y&1) {
        div_i(b,a);
        if(x+y == 3) negate(b,b);
    }
    else if(x+y == 3) negate(b,a);
    else if(&b!=&a) b=a;
}


void add(GG& c, const GG& a, const GG& b) {// c=a+b
    add(c.x, a.x, b.x);
    add(c.y, a.y, b.y);
}

void sub(GG& c, const GG& a, const GG& b) {// c=a-b
    sub(c.x, a.x, b.x);
    sub(c.y, a.y, b.y);
}

void mul(GG& c, const GG& a, const GG& b) {// c=a*b
    ZZ s,t,u,v;
    mul(s, a.x, b.x);
    mul(t, a.y, b.y);
    sub(u, a.y, a.x);
    sub(v, b.x, b.y);
    sub(c.x, s, t);
    mul(c.y, u, v);
    c.y += s;
    c.y += t;
}

void sqr(GG& b, const GG& a) {// b=a*a
    ZZ s,t;
    add(s, a.x, a.y);
    sub(t, a.x, a.y);
    mul(b.y, a.x, a.y);
    mul(b.x, s, t);
    b.y <<= 1;
}

void div(GG& q, const GG& a, const GG& b)
// q = quotient of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2
{
    ZZ n;
    GG c;
    if(IsZero(b.y)) {
        n = b.x;
        c = a;
    }
    else if(IsZero(b.x)) {
        n = b.y;
        div_i(c,a);
    }
    else {
        norm(n,b);
        conj(c,b);
        c *= a;
    }
    c.x <<= 1; c.x += n;
    c.y <<= 1; c.y += n;
    n <<= 1;
    div(q.x, c.x, n);
    div(q.y, c.y, n);
}

void rem(GG& r, const GG& a, const GG& b)
// r = remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2
{
    GG q;
    DivRem(q,r,a,b);
}

void DivRem(GG& q, GG& r, const GG& a, const GG& b)
// q,r = quotient and remainder of a/b such that
//   a = bq + r and |r/b|^2 <= 1/2
{
    if(&a==&q || &a==&r) {
        GG c(a);
        DivRem(q,r,c,b);
        return;
    }
    if(&b==&q) {
        GG c(b);
        DivRem(q,r,a,c);
        return;
    }
    div(q,a,b);
    mul(r,b,q);
    sub(r,a,r);
}

long divide(GG& q, const GG& a, const GG& b)
// if a/b is divisible, set q=a/b and return 1
// else return 0
{
    ZZ n;
    GG c;
    if(IsZero(b.y)) {
        n = b.x;
        c = a;
    }
    else if(IsZero(b.x)) {
        n = b.y;
        div_i(c,a);
    }
    else {
        norm(n,b);
        conj(c,b);
        c *= a;
    }
    if(!divide(c.x, c.x, n) ||
       !divide(q.y, c.y, n)) return 0;
    q.x = c.x;
    return 1;
}

long divide(const GG& a, const GG& b)
// if a/b is divisible, return 1, else return 0
{
    ZZ n;
    GG c;
    if(IsZero(b.y))
        return divide(a.x, b.x) && divide(a.y, b.x);
    if(IsZero(b.x))
        return divide(a.x, b.y) && divide(a.y, b.y);
    norm(n,b);
    conj(c,b);
    c *= a;
    return divide(c.x, n) && divide(c.y, n);
}

void power(GG& b, const GG& a, long n)
// b = a^n; assume n>=0
{
    if(n==0 || IsOne(a)) { set(b); return; }
    if(&b==&a) { GG c(a); power(b,c,n); return; }
    long m(1<<(NumBits(n)-1));
    b=a;
    for(m>>=1; m; m>>=1) {
        sqr(b,b);
        if(n&m) b*=a;
    }
}

long IsUnit(const GG& a) {// return 1 if |a|=1 else 0
    return (IsZero(a.y) && IsOne(a.x) || a.x==-1) ||
           (IsZero(a.x) && IsOne(a.y) || a.y==-1);
}

long IsAssoc(const GG& a, const GG& b)
// return 1 if a = b*i^e for some e (e=0,1,2,3) else 0
{
    if(a==b) return 1;
    GG c;
    mul_i(c,b); if(a==c) return 1;
    mul_i(c,c); if(a==c) return 1;
    mul_i(c,c); if(a==c) return 1;
    return 0;
}

long ProbPrime(const GG& a, long NTRY)
// test if a is gaussian prime, i.e.,
//   |a|^2 is prime and |a|^2==1 (mod 4) or
//   Im(a)==0 and |Re(a)| is prime and |Re(a)|==3 (mod 4) or
//   Re(a)==0 and |Im(a)| is prime and |Im(a)|==3 (mod 4)
// NTRY: number of trials of Miller-Rabin test
{
    ZZ b;
    if(IsZero(a.y)) abs(b, a.x);
    else if(IsZero(a.x)) abs(b, a.y);
    if(!IsZero(b))
        return trunc_long(b,2)==3 && ProbPrime(b, NTRY);
    norm(b,a);
    return b==2 || trunc_long(b,2)==1 && ProbPrime(b, NTRY);
}

void GCD(GG& d, const GG& a, const GG& b)
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
// by euclidean algorithm
{
    GG x(a),y(b),r;
    while(!IsZero(y)) {
        rem(r,x,y);
        x=y;
        y=r;
    }
    FirstQuad(d,x);
}

void XGCD(GG& d, GG& s, GG& t, const GG& a, const GG& b)
// d = greatest common divisor of a and b
//   in first quadrant Re(d)>0 and Im(d)>=0
//   and compute s,t such that d = s*a + t*b
// by extended euclidean algorithm
{
    long c;
    GG x(a),y(b),u,v(1),q,r;
    set(s);
    clear(t);
    while(!IsZero(y)) {
        DivRem(q,r,x,y);
        mul(x,q,u);
        sub(x,s,x);
        s=u;
        u=x;
        mul(x,q,v);
        sub(x,t,x);
        t=v;
        v=x;
        x=y;
        y=r;
    }
    c = FirstQuad(d,x);
    mul_i(s,s,c);
    mul_i(t,t,c);
}

void RandomBnd(GG& a, const ZZ& n)
// a = random gaussian integer such that
//   |Re(a)|<n and |Im(a)|<n
{
    RandomBnd(a.x, n);
    RandomBnd(a.y, n);
    mul_i(a, a, RandomBits_long(2));
}

void GenPrime(GG& p, long l, long f, long err)
// generate random gaussian prime.
// p = gaussian integer such that |p|^2 == q^f
//   where q is random prime integer
// l = bit length of q, so that 2^{l-1} < q < 2^l
// f = 1 or 2 (P == 1 or 3 mod 4, respectively).
//   If f==1, p is imaginary, Re(p)>0, Im(p)>0.
//   If f==2, p = q (real and positive).
// err = bound for error probability < 2^{-err}
// If f<1 or f>2, f is set to 1 or 2, respectively.
{
    if(l<2) Error("l<2 in GenPrime");
    ZZ q;
    if(f<=1) {
        do GenPrime(q,l,err);
        while(trunc_long(q,2)==3);
        if(q==2) { set(p,1,1); return; }
        FactorPrime(p,q);
        if(RandomBits_long(1)) conj(p,p);
    }
    else {
        do GenPrime(p.x, l, err);
        while(trunc_long(p.x, 2)!=3);
        clear(p.y);        
    }
}

std::ostream& operator<<(std::ostream& s, const GG& a) {
    s << '[' << a.x << ' ' << a.y << ']';
    return s;
}
