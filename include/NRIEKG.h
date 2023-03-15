//---------------------------------------------------------------------------
//
// $Id: NRIEKG.h,v 1.5 2004/03/25 20:51:09 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Numerical Recipes functions used by IEKG.
//
// Steve Drasco
//
#ifndef _NRIEKG_H
#define _NRIEKG_H

#include "Globals.h"

Real elle(const Real phi,const Real ak);
Real ellf(const Real phi,const Real ak);
Real ellpi(const Real phi,const Real en,const Real ak);
Real rd(const Real x,const Real y,const Real z);
Real rf(const Real x,const Real y,const Real z);
Real rj(const Real x,const Real y,const Real z,const Real p);

#endif
