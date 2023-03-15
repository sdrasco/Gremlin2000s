//---------------------------------------------------------------------------
//
// $Id: NRCKG.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes for all the Numerical Recipes non-utility functions that
// the Geodesic code uses.
//
// Scott Hughes, 13 November 1998
//
#ifndef _NRCKG_H
#define _NRCKG_H

#include "Globals.h"

Real elle(const Real phi, const Real ak);
Real ellf(const Real phi, const Real ak);
Real ellpi(const Real phi, const Real en, const Real ak);

Real rd(const Real x, const Real y, const Real z);
Real rf(const Real x, const Real y, const Real z);
Real rj(const Real x, const Real y, const Real z, const Real p);

#endif
