from GG import *

# Greatest Common Divisor
a,b,c = GG(2,3),GG(4,5),GG(6,7)
a *= c; b *= c
d = GCD(a,b); print(d)
d,s,t = XGCD(a,b)
print(d == s*a + t*b)

# Factoring integers into Gaussian primes
a,b,c = GenPrime(8),GenPrime(8),GenPrime(8)
a = a**2 * b**3 * c**4; print(a)
f = factor(a); print(f)
b = 1
for k in f: b *= k**f[k]
print(IsAssoc(a,b))

# Biquadratic Residue Symbol
a = GenPrime(20)
p = GenPrime(20); a%=p
s = ResSymb(a,p); print(a,p,s)
print(s == PowerMod(a, (norm(p)-1)//4, p))
if s==1:# solve x^4 == a (mod p)
    x = QrtRootMod(a,p); print(x)
    print(PowerMod(x,4,p) == a)
