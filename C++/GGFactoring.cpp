// uses NTL
//   http://www.shoup.net/ntl

#include "GGFactoring.h"
using namespace NTL;

void cornacchia(ZZ& x, ZZ& y, const ZZ& p)
// given a prime number p where p==1 (mod 4),
// find x,y such that x^2 + y^2 = p and x>y>0
// by cornacchia algorithm
// reference: H. Wada
//   "Prime Factorization by Computers" (in Japanese) p69
{
    ZZ s,r;
    SqrRoot(s,p);
    sub(r,x=p,1);
    SqrRootMod(y,r,p);
    while(x>s) {
        rem(r,x,y);
        x=y;
        y=r;
    }
}

void factor(Vec<Pair<GG, long> >& f, const GG& a)
// f = factorization of a into gaussian primes
// each element of f is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a
// real factors are inserted first in f (if any)
// imaginary factors are appended after real factors
// real factors are positive and sorted in increasing order
// imaginary factors are in first quadrant and sorted by norm
{
    int i,j,k(0),e(0);
    ZZ s,t;
    GG b(a);
    Vec<Pair<ZZ, long> > g,h;
    f.SetLength(0);
    if(IsZero(a) || IsUnit(a)) return;
    GCD(t, real(b), imag(b));
    real(b) /= t;// primitive part
    imag(b) /= t;
    norm(s,b);
    factor(g,s);
    factor(h,t);// content
    for(i=j=0; i<h.length(); i++) {
        if(trunc_long(h[i].a, 2) == 3) {
            f.SetLength(k+1);// real factor
            f[k].a = h[i].a;
            f[k++].b = h[i].b;
        }
        else { if(i>j) h[j] = h[i]; j++; }
    }
    h.SetLength(j);
    i=j=0;
    if(g.length() && g[0].a == 2) e += g[i++].b;
    if(h.length() && h[0].a == 2) e += h[j++].b*2;
    if(e) {// factor 1+i
        f.SetLength(k+1);
        set(f[k].a, 1, 1);
        f[k++].b = e;
    }
    while(i<g.length() || j<h.length()) {// imaginary factor
        if(j==h.length() || i<g.length() && g[i].a < h[j].a) {
            f.SetLength(k+1);
            cornacchia(s, t, g[i].a);
            set(f[k].a, s, t);
            if(!divide(b, f[k].a)) mirror(f[k].a, f[k].a);
            f[k++].b = g[i++].b;
        }
        else {
            f.SetLength(k+2);
            cornacchia(s, t, h[j].a);
            set(f[k].a, s, t);
            mirror(f[k+1].a, f[k].a);
            f[k+1].b = f[k].b = h[j].b;
            if(i<g.length() && g[i].a == h[j].a) {
                if(divide(b, f[k].a)) f[k].b += g[i++].b;
                else f[k+1].b += g[i++].b;
            }
            j++; k+=2;
        }
    }
}

void mul(GG& a, const Vec<Pair<GG, long> >& f)
// a = product of (gaussian integer)^{exponent} in f
// each element of f is a pair of integer and exponent
{
    int i;
    GG b;
    set(a);
    for(i=0; i<f.length(); i++) {
        power(b, f[i].a, f[i].b);
        a *= b;
    }
}