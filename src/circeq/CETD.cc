//---------------------------------------------------------------------------
//
// $Id: CETD.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial Teukolsky driver
// Scott Hughes, 26 July 2003
//
#include <math.h>
#include "Globals.h"
#include "CETD.h"
#include "GKG.h"
#include "CEKG.h"
#include "SWSH.h"
#include "SNT.h"
#include "CEKR.h"
#include "RRGW.h"
#include "NRUtil.h"

CETD::CETD(const int orbitsense, const Real rad, const Real spin,
	   const Real eps, char out[]) :
  proret(orbitsense), r(rad), a(spin), EPSILON(eps)
{
  rrgw = new RRGW;

  sprintf(outname, "%s", out);
  outfile = fopen(outname, "w");
  fclose(outfile);
}

CETD::~CETD()
{
  delete rrgw;
}

void CETD::Driver(const int lmax)
{
  for (l = 2; l <= lmax; l++) {
    //
    // Only do positive m; can get negative m by symmetry.
    //
    for(m = 1; m <= l; m++) {
      DoHarmonic(l, m);
    }
  }
}

void CETD::Driver(const Real EPS_L)
{
  int dE_l_small = 0;
  Real dE_l, dE_l_max = 0;

  l = 2; 
  do { // l-loop
    dE_l = 0.;
    //
    // Only do positive m; can get negative m by symmetry.
    //
    for(m = 1; m <= l; m++) {
      DoHarmonic(l, m);
      dE_l += 2.*(outdata.EdotH + outdata.EdotInf);
    }
    if (dE_l > dE_l_max) dE_l_max = dE_l;
    //
    // The l-loop stops when it's had 3 iterations in a row giving a
    // total energy flux EPS_L smaller than the maximum energy flux
    //
    if (dE_l < EPS_L * dE_l_max)
      dE_l_small++;
    else
      dE_l_small = 0;
    l++;
  } while(dE_l_small <= 2);
}

void CETD::DoHarmonic(const int l, const int m)
{
  CEKG cekg(m, proret, r, a);
  SWSH swsh(l, -2, m, a*cekg.wm);
  SNT snt(m, r, a, cekg.wm, swsh.lambda, EPSILON);
  CEKR cekr(&swsh, &snt, &cekg);
  //
  // Store various quantities in data structure.
  //
  outdata.l = l; outdata.m = m; outdata.k = 0;
  outdata.r = r; outdata.a = a;
  outdata.Lz = cekg.Lz; outdata.E = cekg.E; outdata.Q = 0.;
  outdata.w = cekg.wm; outdata.p = cekg.pm;
  outdata.lambda = swsh.lambda;
  outdata.ZH = cekr.ZedH; outdata.ZInf = cekr.ZedInf;
  rrgw->Flux_Infinity(a, m, swsh.lambda, cekg.wm, cekg.pm, cekr.ZedH,
		      outdata.EdotInf, outdata.LzdotInf);
  rrgw->Flux_Horizon(a, m, swsh.lambda, cekg.wm, cekg.pm, cekr.ZedInf,
		     outdata.EdotH, outdata.LzdotH);
  rrgw->Qdotrdot(r, a, 0.0, cekg.E, cekg.Lz,
		 outdata.EdotInf, outdata.LzdotInf,
		 outdata.QdotInf, outdata.rdotInf);
  rrgw->Qdotrdot(r, a, 0.0, cekg.E, cekg.Lz,
		 outdata.EdotH, outdata.LzdotH,
		 outdata.QdotH, outdata.rdotH);
  outfile = fopen(outname, "a");
  fwrite(&outdata, sizeof(DataHolder), 1, outfile);
  fclose(outfile);
  cout << r << " " << l << " " << m << " " << cekg.wm << " "
       << outdata.EdotH << " " << outdata.EdotInf << endl;
}
