//---------------------------------------------------------------------------
//
// $Id: NRSWSH.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes for all the Numerical Recipes non-utility functions that
// the SWSH code uses.
//
// Scott Hughes, 24 July 1998
//
#ifndef _NRSWSH_H
#define _NRSWSH_H

#include "Globals.h"

Real factrl(const int n);
Real gammln(const Real xx);
Real pythag(const Real a, const Real b);

void indexx(const unsigned long n, Real arr[], unsigned long indx[]);
void tqli(Real d[], Real e[], const int n, Real **z);
void tred2(Real **a, const int n, Real d[], Real e[]);

#endif
