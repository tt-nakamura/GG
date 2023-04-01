#ifndef __GGFactoring_h__
#define __GGFactoring_h__

#include<NTL/pair.h>
#include "GG.h"

void factor(Vec<Pair<ZZ, long> >& f, const ZZ& n);
void factor(Vec<Pair<GG, long> >& f, const GG& a);
void mul(ZZ& a, const Vec<Pair<ZZ, long> >& f);
void mul(GG& a, const Vec<Pair<GG, long> >& f);

#endif // __GGFactoring_h__
