// uses NTL
//   http://www.shoup.net/ntl

#include "GG.h"
using namespace NTL;

std::ostream& operator<<(std::ostream& s, const GG& a) {// print a as [a.x a.y]
    s << '[' << a.x << ' ' << a.y << ']';
    return s;
}

void conv(GG& b, const ZZ& a) { b.x = a; clear(b.y); }// b=a+0i
void conv(GG& b, long a) { b.x = a; clear(b.y); }// b=a+0i

GG& GG::operator=(const ZZ& a) { x=a; clear(y); return *this; }// a+0i
GG& GG::operator=(long a) { x=a; clear(y); return *this; }// a+i0

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

void mul_i(GG& b, const GG& a, long k) {// b = a * i^k
    if((k&=3)==1) mul_i(b,a);
    else if(k==2) negate(b,a);
    else if(k==3) div_i(b,a);
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
// b = a * i^k (k=0,1,2,3) such that Re(b)>0 and Im(b)>=0
// return k; k==0 if a==0
{
    long k(quadrant(a));
    if(k>0) k=4-k; else k=0;
    mul_i(b,a,k);
    return k;
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
// else return 0 (and q is unchanged)
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

long divide2(GG& q, const GG& a)
// if a is divisible by 1+i, set q=a/(1+i) and return 1
// else return 0 and q is unchanged
{
    if(bit(a.x, 0) != bit(a.y, 0)) return 0;
    ZZ x(a.x);
    add(q.x, x, a.y); q.x >>= 1;
    sub(q.y, a.y, x); q.y >>= 1;
    return 1;
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
    return (IsZero(a.y) && (IsOne(a.x) || a.x==-1)) ||
           (IsZero(a.x) && (IsOne(a.y) || a.y==-1));
}

long IsAssoc(const GG& a, const GG& b)
// return 1 if a = b*i^k for some k (k=0,1,2,3) else 0
{
    if(a==b) return 1;
    GG c;
    mul_i(c,b); if(a==c) return 1;
    mul_i(c,c); if(a==c) return 1;
    mul_i(c,c); if(a==c) return 1;
    return 0;
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

// a = random gaussian integer in square [-L,L]^2
void RandomBits(GG& a, long l)// 0 <= L < 2^l
{ RandomBits(a.x, l); RandomBits(a.y, l); mul_i(a, a, RandomBits_long(2)); }

void RandomLen(GG& a, long l)// 2^{l-1} <= L < 2^l
{ RandomLen(a.x, l); RandomLen(a.y, l); mul_i(a, a, RandomBits_long(2)); }

void RandomBnd(GG& a, const ZZ& n)// 0 <= L < n
{ RandomBnd(a.x, n); RandomBnd(a.y, n); mul_i(a, a, RandomBits_long(2)); }

void PowerMod(GG& b, const GG& a, long n, const GG& m)
// b = a^n mod m; assume n>=0 and |a| < |m|
{
    if(n==0 || IsOne(a)) { set(b); return; }
    if(&b==&a) { GG c(a); PowerMod(b,c,n,m); return; }
    long k(1<<(NumBits(n)-1));
    b=a;
    for(k>>=1; k; k>>=1) {
        sqr(b,b); b%=m;
        if(n&k) { b*=a; b%=m; }
    }
}

void PowerMod(GG& b, const GG& a, const ZZ& n, const GG& m)
// b = a^n mod m; assume n>=0 and |a| < |m|
{
    if(IsZero(n) || IsOne(a)) { set(b); return; }
    if(&b==&a) { GG c(a); PowerMod(b,c,n,m); return; }
    b=a;
    for(long k=NumBits(n)-2; k>=0; k--) {
        sqr(b,b); b%=m;
        if(bit(n,k)) { b*=a; b%=m; }
    }
}

long ProbPrime(const GG& a, long NTRY)
// test if a is gaussian prime, i.e., return 1 if either
//   |a| is prime and |a|==3 (mod 4) or
//   |a|^2 is prime and |a|^2==1 (mod 4) or
//   |a|^2 == 2
// else return 0
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

void GenPrime(GG& p, long l, long f, long err)
// generate random gaussian prime p.
// f must be 1 or 2; l must be l>=2
// if f==1, |p|^2=q and q==1 (mod 4)
// if f==2, p=-q+0i and q==3 (mod 4)
// where q is random prime and 2^{l-1} <= q < 2^l
// probability of error is less than 2^-err
// p=x+iy is primary (x odd, y even, x+y==1 (mod 4))
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
        do GenPrime(q,l,err);
        while(trunc_long(q,2)!=3);
        negate(p.x, q);
        clear(p.y);        
    }
}

long primary(GG& b, const GG& a)
// b = unit * a = x+iy such that
//   x==1 and y==0 or x==3 and y==2 (mod 4)
// return k such that a = i^k * b (k=0,1,2,3)
// Assume |a|^2 is odd
{
    long x(a.x%4), y(a.y%4);
    if(y&1) {
        div_i(b,a);
        if(x+y != 3) return 1;
        negate(b,b);
        return 3;
    }
    if(x+y == 3) {
        negate(b,a);
        return 2;
    }
    if(&b!=&a) b=a;
    return 0;
}

void ResSymb(GG& s, const GG& a, const GG& b)
// s = biquadratic residue symbol (a/b)_4 = 0,1,i,-1,-i
// Assume |a| < |b| and |b|^2 is odd
// Assume b is primary, but may not be prime
// reference: K. S. Williams
//   "On the Supplement to the Law of Biquadratic Reciprocity"
//   Proceedings of the American Mathematical Society 59 (1976) 19
{
    long j(0),k,m,n;
    GG u(a),v(b),w;
    while(!IsZero(u)) {
        m = trunc_long(v.x, 4);// m odd
        n = trunc_long(v.y, 4);// n even
        if(sign(v.x) < 0) m = -m;
        if(sign(v.y) < 0) n = -n;
        k = (m-n)>>2;
        if(n &= 2) k--;
        m >>= 1; m &= 3; k &= 3;

        while(divide2(u,u)) j += k;// supplementary law for 1+i
        if(k = primary(u,u)) j -= k*m;// supplementary law for units
        if(n && bit(u.y,1)) j += 2;// reciprocity
        j &= 3;
        
        rem(w,v,u);
        v = u;
        u = w;
    }
    if(!IsUnit(v)) clear(s);
    else if(j==0) set(s);
    else if(j==1) set(s,0,1);
    else if(j==2) conv(s,-1);
    else set(s,0,-1);
}