//---------------------------------------------------------------------------
//
// $Id: CEKG.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial Kerr geodesics.
// Scott Hughes, 11 January 1999.
//
#include <math.h>
#include "Globals.h"
#include "CEKG.h"
#include "NRUtil.h"

CEKG::CEKG(const int emm, const int orbitsense, const Real rad,
	   const Real spin) : m(emm), r(rad), a(spin)
{
  if (orbitsense == PROGRADE) {
    E = Kerr::Eeqpro(r, a);
    Lz = Kerr::Lzeqpro(r, a);
    Om_phi = Kerr::Omega_phi_eqpro(r, a);
  } else {
    E = Kerr::Eeqret(r, a);
    Lz = Kerr::Lzeqret(r, a);
    Om_phi = Kerr::Omega_phi_eqret(r, a);
  }
  wm = m*Om_phi;
  pm = wm - m*a/(2.*Kerr::rplus(a));
}
