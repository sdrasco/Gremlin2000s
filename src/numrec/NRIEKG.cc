//
// Prototypes for all the Numerical Recipes non-utility functions that
// the Geodesic code uses. Taken from NRCKG.h
//
// Steve Drasco, 30 June 2003
//
#ifndef _NRIEKG_H
#define _NRIEKG_H

#include "NRUtil.h"

Real elle(const Real phi, const Real ak);
Real ellf(const Real phi, const Real ak);
Real ellpi(const Real phi, const Real en, const Real ak);

Real rd(const Real x, const Real y, const Real z);
Real rf(const Real x, const Real y, const Real z);
Real rj(const Real x, const Real y, const Real z, const Real p);

#endif
