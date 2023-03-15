//---------------------------------------------------------------------------
//
// $Id: CKR.h,v 1.2 2004/09/04 19:28:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Class CKR (Circular Kerr Radiation).
//
// Scott Hughes, culled from old class Rad_Circ, 9 January 1999.
//
#ifndef _CKR_H
#define _CKR_H

#include "Globals.h"
#include "SWSH.h"
#include "GKG.h"
#include "CKG.h"
#include "SNT.h"

class CKR {
public:
  CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon);

  Complex ZedH;
  Complex ZedInf;

private:
  Complex rho_func(const int n, const Real z);
  Complex Cnn(const int seg, const int n, const Real z);
  Complex Cnmbar(const int seg, const int n, const Real z);
  Complex Cmbarmbar(const int seg, const int n, const Real z);
  Complex Ann0(const int seg, const int n, const Real z);
  Complex Anmbar0(const int seg, const int n, const Real z);
  Complex Ambarmbar0(const int seg, const int n, const Real z);
  Complex Anmbar1(const int seg, const int n, const Real z);
  Complex Ambarmbar1(const int seg, const int n, const Real z);
  Complex Ambarmbar2(const int seg, const int n, const Real z);

  Real Psi_mk(const Real chi);

  // The above methods are used to put together ...
  Complex dZedH(const int seg, const int n, const Real chi);
  Complex dZedInf(const int seg, const int n, const Real chi);

  // Wrapper for the dZeds.
  Complex IntegrandH(const Real chi), IntegrandInf(const Real chi);

  // Wrapper for all the C-C quadrature code used to compute ZedH and
  // ZedInf.
  void ComputeZeds(Complex & ZedH, Complex & ZedInf, const Real epsilon);

  // Code for the cosine fft.
  void cosft1(Complex y[], int n);
  void realft(Complex data[], unsigned long n, int isign);
  void four1(Complex data[], unsigned long nn, int isign);

  // And for the Richardson extrapolation.
  void ratint(const Real xa[], const Complex ya[], const int n,
	      const Real x, Complex *y, Complex *dy);

  // Other stuff we need in lots of places.
  int l, m, k;
  SWSH *swsh;
  SNT *snt;
  GKG gkg;
  CKG *ckg;
  Real r, a, E, Lz, Q;
  Real wmk, pmk;
  Real Delta, dr_Delta;
  Real SN_Kay, dr_SN_Kay;
};

#endif
