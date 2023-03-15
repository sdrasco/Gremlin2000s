//---------------------------------------------------------------------------
//
// $Id: CEKR.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Class CEKR (Circular, Equatorial Kerr Radiation).
//
// Scott Hughes, 11 January 1999.
//
#ifndef _CEKR_H
#define _CEKR_H

#include "Globals.h"
#include "SWSH.h"
#include "GKG.h"
#include "CEKG.h"
#include "SNT.h"

class CEKR {
public:
  CEKR(SWSH *swsh_in, SNT *snt_in, CEKG *cekg_in);

  Complex ZedH;
  Complex ZedInf;

private:
  Complex rho_func();
  Complex Cnn();
  Complex Cnmbar();
  Complex Cmbarmbar();
  Complex Ann0();
  Complex Anmbar0();
  Complex Ambarmbar0();
  Complex Anmbar1();
  Complex Ambarmbar1();
  Complex Ambarmbar2();

  // Other stuff we need in lots of places.
  int l, m;
  SWSH *swsh;
  SNT *snt;
  GKG gkg;
  CEKG *cekg;
  Real r, a, E, Lz;
  Real wm, pm;
  Real Delta, dr_Delta, ddr_Delta;
  Real SN_Kay, dr_SN_Kay;
  Real r_plus;

  // The above methods are used to put together ...
  Complex Zed(const Complex R,
	      const Complex dr_R,
	      const Complex ddr_R);

  //
  // The following functions develop Zed as published in
  // Poisson 93, for checking my routines in the limit of
  // a = 0.
  //
  Real zerob(), neg1b(), neg2b();
  Complex EricZed(const Complex R,
		  const Complex dr_R,
		  const Complex ddr_R);
};

#endif
