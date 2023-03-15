//---------------------------------------------------------------------------
//
// $Id: NREllf.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real ellf(const Real phi, const Real ak)
{
  Real rf(const Real x, const Real y, const Real z);
  const Real s = sin(phi);

  return s*rf(SQR(cos(phi)), (1.0-s*ak)*(1.0+s*ak), 1.0);
}
