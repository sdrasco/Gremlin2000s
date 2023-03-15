//---------------------------------------------------------------------------
//
// $Id: CEKG.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial Kerr geodesics.
// Scott Hughes, 11 January 1999.
//
#ifndef _CEKG_H
#define _CEKG_H

#define PROGRADE 1
#define RETROGRADE 0

#include <math.h>
#include "Globals.h"

class CEKG {
public:
  CEKG(const int emm, const int orbitsense, const Real rad, const Real spin);

  int m;
  Real r, a;

  Real E, Lz, Om_phi;

  Real wm, pm;
};
#endif
