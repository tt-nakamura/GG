class GG:
    """ Gaussian integer x+iy """
    def __init__(a, x=0, y=0):
        a.x = int(x)
        a.y = int(y)

    def __repr__(a):
        if a.y == 0: return str(a.x)
        elif a.x == 0: return str(a.y) + 'j'
        elif a.y < 0:
            return '(' + str(a.x) + str(a.y) + 'j)'
        else:
            return '(' + str(a.x) + '+' + str(a.y) + 'j)'

    def __hash__(a):
        return a.x ^ a.y

    def __neg__(a):
        return GG(-a.x, -a.y)

    def __add__(a,b):
        if isinstance(b, GG):
            return GG(a.x + b.x, a.y + b.y)
        else:
            return GG(a.x + b, a.y)

    def __sub__(a,b):
        if isinstance(b, GG):
            return GG(a.x - b.x, a.y - b.y)
        else:
            return GG(a.x - b, a.y)

    def __mul__(a,b):
        if isinstance(b, GG):
            s = a.x*b.x
            t = a.y*b.y
            u = (a.x - a.y)*(b.y - b.x)
            return GG(s-t, s+t+u)
        else:
            return GG(b*a.x, b*a.y)

    def __truediv__(a,b):
        """ a = bq + r, norm(r) <= norm(b)/2 """
        if isinstance(b, GG):
            a *= conj(b)
            b = norm(b)

        x = ((a.x<<1) + b)//(b<<1)
        y = ((a.y<<1) + b)//(b<<1)
        return GG(x,y)

    def __mod__(a,b):
        return a - a/b*b

    def __divmod__(a,b):
        q = a/b
        return q, a - b*q

    def __pow__(a,n):
        """ assume n>=0 """
        if n==0: return GG(1)
        k = (1<<(n.bit_length()-1))>>1
        b = a
        while k:
            b *= b
            if n&k: b *= a
            k>>=1
        return b

    def __eq__(a,b):
        if isinstance(b, GG):
            return a.x == b.x and a.y == b.y
        else:
            return a.x == b and a.y == 0

    def __radd__(a,b):
        return GG(b + a.x, a.y)

    def __rsub__(a,b):
        return GG(b - a.x, -a.y)

    def __rmul__(a,b):
        return GG(b*a.x, b*a.y)

    def __rtruediv__(a,b):
        return GG(b).__truediv__(a)

    def __rmod__(a,b):
        return GG(b).__mod__(a)

    def __rdivmod__(a,b):
        return GG(b).__divmod__(a)

    def __complex__(a):
        return complex(a.x, a.y)

    def isUnit(a):
        return (a.y==0 and abs(a.x)==1) or\
               (a.x==0 and abs(a.y)==1)

def real(a): return a.x
def imag(a): return a.y
def norm(a): return a.x**2 + a.y**2
def conj(a): return GG(a.x, -a.y)
def mul_i(a): return GG(-a.y, a.x) # a*i
def div_i(a): return GG(a.y, -a.x) # a/i

def IsAssoc(a, b, exponent=False):
    """ a,b: GG, return int
    test if a == b*i^k for some k (k=0,1,2,3)
    if exponent is False, return True or False
    else if a is associate of b, return k
    else return -1
    """
    if   a == b: k=0
    elif a == -b: k=2
    elif a == mul_i(b): k=1
    elif a == div_i(b): k=3
    else: k=-1
    if exponent: return k
    else: return k>=0

def mul_ipow(a,k): # a*i^k
    """ a: GG, k: int, return GG """
    k &= 3
    if   k==2: return -a
    elif k==1: return mul_i(a)
    elif k==3: return div_i(a)
    else: return a

def quadrant(a):
    """ a: GG, return int
    return -1 if a==0,
            0 if Re(a)>0 and Im(a)>=0,
            1 if Re(a)<=0 and Im(a)>0,
            2 if Re(a)<0 and Im(a)<=0,
            3 if Re(a)>=0 and Im(a)<0
    """
    if a.x==0:
        if a.y==0: return -1
        elif a.y>0: return 1
        else: return 3
    elif a.x>0:
        if a.y<0: return 3
        else: return 0
    elif a.y>0: return 1
    else: return 2

def FirstQuad(a, exponent=False):
    """ a: GG, return GG
    return b = a*i^k for some k (k=0,1,2,3)
    in first quadrant Re(b)>0 and Im(b)>=0
    if exponent is True, return b and k
    if a is zero, then b=0 and k=0
    """
    k = quadrant(a)
    k = 4-k if k>0 else 0
    b = mul_ipow(a,k)
    if exponent: return b,k
    else: return b

def GCD(a,b):
    """ a,b: GG, return GG
    d = greatest common divisor of a and b
      in first quadrant Re(d)>0 and Im(d)>=0
      by Euclidean algorithm
    return d
    if a and b are both zero, then d=0
    """
    while b!=0: a,b = b,a%b
    return FirstQuad(a)

def XGCD(a,b):
    """ a,b: GG, return GG*3
    d = greatest common divisor of a and b
      in first quadrant Re(d)>0 and Im(d)>=0
      by extended Euclidean algorithm
    return d,s,t such that d = s*a + t*b
    if a and b are both zero, then d,s,t=0,1,0
    """
    s,t = GG(1),GG(0)
    u,v = GG(0),GG(1)
    while b!=0:
        q,r = divmod(a,b)
        s,u = u, s-u*q
        t,v = v, t-v*q
        a,b = b,r
    d,e = FirstQuad(a, True)
    s = mul_ipow(s,e)
    t = mul_ipow(t,e)
    return d,s,t

def primary(a, exponent=False):
    """ a: GG, return GG
    return b=x+iy such that a = b*i^k (k=0,1,2,3)
      and x is odd and y is even and x+y==1 (mod 4).
    assume norm(a) is odd
    if exponent is True, return b and k
    """
    x,y = a.x&3, a.y&3
    if x+y==3: a,k = -a,2
    else: k = 0
    if y&1: a,k = div_i(a), k+1
    if exponent: return a,k
    else: return a

def ResSymb(a,b):
    """a,b: GG, return GG
    return biquadratic residue symbol (a/b)_4 = 0,1,i,-1,-i
    Assume norm(a) < norm(b) and norm(b) is odd
    Assume b is primary, but may not be prime
    reference: K. S. Williams
      "On the Supplement to the Law of Biquadratic Reciprocity"
      Proceedings of the American Mathematical Society 59 (1976) 19
    """
    j,lam = 0, GG(1,1)
    while a!=0:
        m = (((b.x - b.y)>>2) - ((b.y>>1)&1))&3
        n = (b.x>>1)&3
        while True:
            q,r = divmod(a,lam)
            if r!=0: break
            a,j = q, j+m

        a,k = primary(a, True)
        if a.y&2 and b.y&2: j += 2
        a,b,j = b%a, a, (j-k*n)&3

    if not b.isUnit(): return 0
    elif j==0: return GG(1)
    elif j==1: return GG(0,1)
    elif j==2: return GG(-1)
    else: return GG(0,-1)

def PowerMod(a,n,m): # a^n mod m
    """ a,m: GG, n: int, return GG
    assume n>=0 """
    if n==0: return GG(1)
    k = (1<<(n.bit_length()-1))>>1
    b = a
    while k:
        b = b*b%m
        if n&k: b = b*a%m
        k>>=1
    return b

#########################################################
import sympy as sp
from sympy.abc import x
from math import gcd
from random import randrange

def SqrRoot(n): # floor(sqrt(n))
    """ n: int, return int """
    x,y = n+1,n
    while x>y: x,y = y, (y + n//y)//2
    return x

def FactorPrime(p):
    """ p: int, return GG
    given prime number p,
    find x,y such that x*x + y*y = p and
      x odd, y even, x+y==1 (mod 4)
    return primary(x+yi)
    Assume p==1 (mod 4)
    referece: H. Wada
      "Prime Factorization by Computers" (in Japanese) p69
    """
    x,y = p, sp.sqrt_mod(p-1, p)
    s = SqrRoot(p)
    while x>s: x,y = y,x%y
    return primary(GG(x,y))

def factor(a):
    """ a: GG, return dict{GG,int}
    factorize a into gaussian primes
    return dictionary of (prime, exponent) pair
    such that product of p**e is associate of a
    real factors are positive and
    imaginary factors are primary
    """
    if not isinstance(a,GG): a = GG(a)
    f = {}
    if a==0 or a.isUnit(): return f
    d = gcd(real(a), imag(a))
    if d>1: a /= d
    g = sp.factorint(norm(a))
    h = sp.factorint(d)

    k = g.pop(2,0) + 2*h.pop(2,0)
    if k: f[GG(1,1)] = k
    for k in h:
        if k&3 == 3: f[k] = h[k]
        else:
            p = FactorPrime(k)
            q = conj(p)
            f[p],f[q] = h[k],h[k]
            if k in g:
                if a%p == 0: f[p] += g.pop(k)
                else:        f[q] += g.pop(k)
    for k in g:
        p = FactorPrime(k)
        if a%p != 0: p = conj(p)
        f[p] = g[k]

    return f

def IsPrime(a):
    """ a: GG, return bool
    test if a is gaussian prime """
    if not isinstance(a,GG): a = GG(a)
    if   imag(a)==0: b = abs(real(a))
    elif real(a)==0: b = abs(imag(a))
    else: b = -1
    if b>=0: return b&3 == 3 and sp.isprime(b)
    b = norm(a)
    return b==2 or (b&3 == 1 and sp.isprime(b))

def GenPrime(l, f=1):
    """ l,f: int, return GG
    generate random gaussian prime p
    Assume f=1 or 2; Assume l>=2
    if f=1, |p|^2=q and q==1 (mod 4)
    if f=2, p=-q+0i and q==3 (mod 4)
    where q is random prime and 2^{l-1} <= q < 2^l
    p=x+iy is primary (x odd, y even, x+y==1 (mod 4))
    """
    while True:
        p = sp.randprime(1<<(l-1), 1<<l)
        if f <= 1:
            if p == 2: return GG(1,1)
            if p&3 == 1:
                p = FactorPrime(p)
                if randrange(2): return p
                else: return conj(p)

        elif p&3 == 3: return GG(-p)

def InvMod(a,p):
    """ a,p: int, return int
    return x such that ax==1 (mod p)
    Assume a and p are relatively prime
    """
    s,u = 1,0
    while p:
        q,r = divmod(a,p)
        a,p,s,u = p,r,u,s-q*u
    return s

def QrtRootMod(a,p):
    """ a,p: GG, return GG
    solve x^4 == a (mod p) and return x
    Assume norm(p) is prime and norm(p)==1 (mod 4)
    Assume (a/p)_4 == 1
    """
    n = norm(p)
    b = p.x * InvMod(p.y, n) % n
    b = (a.x - b * a.y % n) % n
    F = sp.factor_list(x**4 - b, modulus=n)
    return int(-F[1][0][0].coeff(x,0))%p
