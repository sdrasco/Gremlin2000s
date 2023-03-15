//---------------------------------------------------------------------------
//
// $Id: CETD.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial Teukolsky driver
// Scott Hughes, 26 July 2003
//
#ifndef _CETD_H
#define _CETD_H

#include <math.h>
#include "Globals.h"
#include "CEKG.h"
#include "SWSH.h"
#include "SNT.h"
#include "CEKR.h"
#include "RRGW.h"
#include "NRUtil.h"

class CETD {
public:
  CETD(const int orbitsense, const Real rad, const Real spin,
       const Real eps, char out[]);
  ~CETD();

  void Driver(const int lmax);
  void Driver(const Real EPS_L);

private:
  void DoHarmonic(const int l, const int m);

  DataHolder outdata;
  char outname[50];
  FILE *outfile;

  RRGW *rrgw;
  int proret;
  int l, m;
  Real r, a, EPSILON;
};
#endif
