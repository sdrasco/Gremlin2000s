//---------------------------------------------------------------------------
//
// $Id: CTD.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Circular Teukolsky driver
// Scott Hughes, 28 July 2003
//
#ifndef _CTD_H
#define _CTD_H

#include <math.h>
#include "Globals.h"
#include "CEKG.h"
#include "CKG.h"
#include "SWSH.h"
#include "SNT.h"
#include "CEKR.h"
#include "CKR.h"
#include "RRGW.h"
#include "NRUtil.h"

#define EPSILON_CKR 1.e-8

class CTD {
public:
  CTD(const Real rad, const Real Lz_or_E_or_cosiota_in,
      const int Lz_or_E_or_cosiota_flag_in, const Real spin,
      const Real eps, char out[]);
  ~CTD();

  void Driver(const Real EPS_L, const Real EPS_K);

private:
  void DoHarmonic(const int l, const int m, const int k);

  DataHolder outdata;
  char outname[50];
  FILE *outfile;

  RRGW *rrgw;
  int l, m, k;
  Real r, Lz_or_E_or_cosiota;
  int Lz_or_E_or_cosiota_flag;
  Real a, EPSILON;
};
#endif
