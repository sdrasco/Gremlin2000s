//---------------------------------------------------------------------------
//
// $Id: CKG.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular Kerr geodesics.
//
// Scott Hughes, pulled from old header 'geodesic.h',
// 10 January 1999.
//
// 26 September 1999: added a flag, Lz_or_E_flag, and changed
// orbitL to Lz_or_E.  This is so that the same class can be used
// to cover the nearly horizon skimming orbits as all other
// orbits.
//
#ifndef _CKG_H
#define _CKG_H

#include <math.h>
#include "Globals.h"
#include "NRCKG.h"

#define USING_cosiota 2
#define USING_Lz 1
#define USING_E 0

class CKG {
public:
  CKG(const int emm, const int kay, const Real Lz_or_E_or_cosiota,
      int Lz_or_E_or_cosiota_flag, const Real rad, const Real spin);

  // These utility constants are used in many of the functions
  // below, and are set by the constructor.
  int m, k;
  Real Lz, E, r, a, Delta;
  Real alpha_wilkins, beta_wilkins;
  Real zedplus, zedminus;
  Real betazedplus;
  Real k_zmzp;   // Parameter k in elliptical integrals.
  Real gamma, delta;
  Real t_pi2, phi_pi2;

  // The following functions and quantities are constant for
  // fixed Lz, r, a.  The constants are set by the constructor.
  Real thetamax, thetamin;
  Real Efunc(), Lzfunc();
  Real Qfunc(), Q;
  Real cosiota;
  Real Period_Theta(), T_th;
  Real Phi(), Phi_range;
  Real Omega_Theta(), Om_th;
  Real Omega_Phi(), Om_phi;

  // Frequency is built from the above methods.
  Real omegamk(), wmk, pmk;

  // The following quantities depend on z.
  //  Real dtdz(const Real z);  // The function {\cal Z} in my notes.
  //
  // Note dtdz has been dropped in favor of dtdchi.
  Real dtdchi(const Real z);
  Real t_0(const Real chi);
  Real t(const Real chi);
  Real phi_0(const Real chi);
  Real phi(const Real chi);
};
#endif
