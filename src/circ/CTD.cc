//---------------------------------------------------------------------------
//
// $Id: CTD.cc,v 1.2 2004/09/04 19:28:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial Teukolsky driver
// Scott Hughes, 26 July 2003
//
#include <math.h>
#include "Globals.h"
#include "CTD.h"
#include "GKG.h"
#include "CEKG.h"
#include "SWSH.h"
#include "SNT.h"
#include "CEKR.h"
#include "RRGW.h"
#include "NRUtil.h"

CTD::CTD(const Real rad, const Real Lz_or_E_or_cosiota_in,
	 int Lz_or_E_or_cosiota_flag_in, const Real spin, const Real eps,
	 char out[]) :
  r(rad), Lz_or_E_or_cosiota(Lz_or_E_or_cosiota_in),
  Lz_or_E_or_cosiota_flag(Lz_or_E_or_cosiota_flag_in),
  a(spin), EPSILON(eps)
{
  rrgw = new RRGW;

  sprintf(outname, "%s", out);
  outfile = fopen(outname, "w");
  fclose(outfile);
}

CTD::~CTD()
{
  delete rrgw;
}

void CTD::Driver(const Real EPS_L, const Real EPS_K)
{
  int dE_l_small = 0;
  Real dE_l, dE_l_max = 0;
  int dE_lmk_small = 0;
  Real dE_lmk, dE_lmk_max = 0., prev_dE_lmk = 0.;
  //
  l = 2; 
  do { // l-loop
    dE_l = 0.;
    //
    // Do all m, k >= 0; can get rest by symmetry
    //
    for(m = -l; m <= l; m++) {
      if (m > 0) k = 0;
      else k = 1;
      dE_lmk_max = 0.;
      dE_lmk_small = 0;
      do { // k-loop
	dE_lmk = 0.;
	DoHarmonic(l, m, k);
	dE_lmk = 2.*(outdata.EdotH + outdata.EdotInf);
	if (dE_lmk > dE_lmk_max) dE_lmk_max = dE_lmk;
	//
	// Only allow the possibility of terminating if wmk > 0 ...
	// otherwise, we could permaturely kill the run and miss a
	// subpeak of the distribution.
	//
	// The k-loop stops when it's had 3 iterations in a row giving
	// a total energy flux EPS_K smaller than the maximum energy
	// flux in this k-sequence
	//
	if (outdata.w > 0) {
	  if (fabs(dE_lmk) <= EPS_K * dE_lmk_max &&
	      fabs(dE_lmk) <= fabs(prev_dE_lmk)) {
	    dE_lmk_small++;
	  } else {
	    dE_lmk_small = 0;
	  }
	}
	k++;
	dE_l += dE_lmk;
	prev_dE_lmk = dE_lmk;
      } while(dE_lmk_small <= 2);
    } // m-loop
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

void CTD::DoHarmonic(const int l, const int m, const int k)
{
  CKG ckg(m, k, Lz_or_E_or_cosiota, Lz_or_E_or_cosiota_flag, r, a);
  SWSH swsh(l, -2, m, a*ckg.wmk);
  //
  // If wmk is very very very small, the next step will take a month
  // and the result will be tiny.  However, we don't want to simply
  // output nothing --- we should at least put a blank space with the
  // frequency of the mode in the output file.
  //
  if (fabs(ckg.wmk*ckg.T_th) > EPSILON) {
    SNT snt(m, r, a, ckg.wmk, swsh.lambda, Min(EPSILON, 1e-6));
    //
    // Pass them to the radiation calculator.
    //
    CKR ckr(&swsh, &snt, &ckg, EPSILON_CKR);
    //
    // Store various quantities in data structure.
    //
    outdata.l = l; outdata.m = m; outdata.k = k;
    outdata.r = r; outdata.a = a;
    outdata.Lz = ckg.Lz; outdata.E = ckg.E;
    outdata.Q = ckg.Q; outdata.w = ckg.wmk;
    outdata.p = ckg.pmk; outdata.lambda = swsh.lambda;
    //
    outdata.ZH = ckr.ZedH;
    outdata.ZInf = ckr.ZedInf;
    rrgw->Flux_Infinity(a, m, outdata.lambda, outdata.w,
			outdata.p, outdata.ZH,
			outdata.EdotInf, outdata.LzdotInf);
    rrgw->Flux_Horizon(a, m, outdata.lambda, outdata.w,
		       outdata.p, outdata.ZInf,
		       outdata.EdotH, outdata.LzdotH);
    rrgw->Qdotrdot(r, a, outdata.Q, outdata.E, outdata.Lz,
		   outdata.EdotInf, outdata.LzdotInf,
		   outdata.QdotInf, outdata.rdotInf);
    rrgw->Qdotrdot(r, a, outdata.Q, outdata.E, outdata.Lz,
		   outdata.EdotH, outdata.LzdotH,
		   outdata.QdotH, outdata.rdotH);
  } else {
    outdata.l = l; outdata.m = m; outdata.k = k;
    outdata.r = r; outdata.a = a;
    outdata.Lz = ckg.Lz; outdata.E = ckg.E;
    outdata.Q = ckg.Q; outdata.w = ckg.wmk;
    outdata.p = ckg.pmk; outdata.lambda = swsh.lambda;
    //
    outdata.ZH = 0.0;
    outdata.ZInf = 0.0;
    outdata.EdotInf = 0.0; outdata.LzdotInf = 0.0;
    outdata.EdotH = 0.0; outdata.LzdotH = 0.0;
    outdata.QdotInf = 0.0; outdata.rdotInf = 0.0;
    outdata.QdotH = 0.0; outdata.rdotH = 0.0;
  }
  outfile = fopen(outname, "a");
  fwrite(&outdata, sizeof(DataHolder), 1, outfile);
  fclose(outfile);
  //cout << r << " " << l << " " << m << " " << k << " " << ckg.wmk << " "
  //     << outdata.EdotH << " " << outdata.EdotInf << " " << endl;
  cout.precision(16);
  cout << l << " " << m << " " << k << " " << ckg.wmk << " "
       << outdata.EdotH << " " << outdata.EdotInf << " " 
       << outdata.LzdotH << " " << outdata.LzdotInf << endl;

}
