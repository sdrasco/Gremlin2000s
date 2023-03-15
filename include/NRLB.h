//---------------------------------------------------------------------------
//
// $Id: NRLB.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes for the numerical recipes functions needed by the LB code.
//
// Scott Hughes, 19 September 1999
//
#ifndef _NRLB_H
#define _NRLB_H

#include "Globals.h"

void lubksb(Real **a, int n, int *indx, Real b[]);
void ludcmp(Real **a, int n, int *indx, Real *d);

#endif
