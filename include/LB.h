//---------------------------------------------------------------------------
//
// $Id: LB.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes and class structure for the LB code.
//
// Scott Hughes, 3 September 1999
//
#ifndef _LB_H
#define _LB_H

#include "Globals.h"
#include "GKG.h"

class LB {
public:
  LB(const Real arr, const Real ayy);
  ~LB();

  Real E, Lz, Q, cosiota;

private:
  void func_and_jac(Real *x, int n, Real *fvec, Real **fjac);
  void mnewt(const int ntrial, Real x[], const int n, const Real tolx,
	     const Real tolf);

  Real r, a;
  GKG *gkg;
};

#endif
