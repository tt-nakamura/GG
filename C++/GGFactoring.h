// uses NTL
//   http://www.shoup.net/ntl

#ifndef __GGFactoring_h__
#define __GGFactoring_h__

#include<NTL/pair.h>
#include "GG.h"

void factor(NTL::Vec<NTL::Pair<NTL::ZZ, long> >& f, const NTL::ZZ& n);
// f = prime factorization of |n|
//   vector of (prime, exponent) pair
//   in increasing order of primes

void factor(NTL::Vec<NTL::Pair<GG, long> >& f, const GG& a);
// f = factorization of a into gaussian primes
// each element of f is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a
// real factors are inserted first in f (if any)
// imaginary factors are appended after real factors
// real factors are positive and sorted in increasing order
// imaginary factors are in first quadrant and sorted by norm

void mul(GG& a, const NTL::Vec<NTL::Pair<GG, long> >& f);
// a = product of (gaussian integer)^{exponent} in f
// each element of f is a pair of integer and exponent

#endif // __GGFactoring_h__
