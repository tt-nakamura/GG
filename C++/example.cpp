#include "GGFactoring.h"
#include<iostream>
using namespace NTL;

main() {
    GG a(2,3),b(4,5),c(6,7);
    GG d,s,t;

    // Greatest Common Divisor
    a *= c; b *= c;
    GCD(d,a,b); std::cout << d << std::endl;
    XGCD(d,s,t,a,b);
    s *= a; t *= b; s += t;
    std::cout << (d==s) << std::endl;

    // Factoring integers into Gaussian primes
    Vec<Pair<GG,long> > f;
    GenPrime(a,8);
    GenPrime(b,8);
    GenPrime(c,8);
    power(a,a,2);
    power(b,b,3); a *= b;
    power(c,c,4); a *= c;
    factor(f,a); std::cout << f << std::endl;
    mul(b,f); std::cout << IsAssoc(a,b) << std::endl;

    // Biquadratic Residue Symbol
    GG p;
    ZZ n;
    GenPrime(p,20);
    GenPrime(a,20); a%=p;
    ResSymb(c,a,p); std::cout << c << std::endl;
    norm(n,p); n--; n >>= 2;
    PowerMod(b,a,n,p);
    std::cout << (b==c) << std::endl;
    if(IsOne(c)) {// solve x^4 == a (mod p)
        QrtRootMod(b,a,p); std::cout << b << std::endl;
        PowerMod(b,b,4,p);
        std::cout << (a==b) << std::endl;
    }
}