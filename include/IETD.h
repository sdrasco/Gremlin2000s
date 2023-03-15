//
// Inclined-EccentricTeukolsky driver
// Steve Drasco, 3 January 2006
//
#ifndef _IETD_H
#define _IETD_H

#include <math.h>
#include "Globals.h"
#include "CEKG.h"
#include "IEKG.h"
#include "CKG.h"
#include "SWSH.h"
#include "SNT.h"
#include "CEKR.h"
#include "CKR.h"
#include "RRGW.h"
#include "NRUtil.h"

class IETD {
public:
  IETD(IEKG *iekg_in, Real eps_in);
  ~IETD();

  // input
  IEKG  *iekg; // Inclined Eccentric Kerr Geodesic object

  // constants set by the constructor
  Real dEH;  // <dE/dt>_H       (M/mu)^2
  Real dEI;  // <dE/dt>_\infty  (M/mu)^2
  Real dLH;  // <dL/dt>_H       M/mu^2
  Real dLI;  // <dL/dt>_\infty  M/mu^2
  Real dQ;   // <dQ/dt>         / (M mu^2)
  Real eps; 

private:

  int MAXSTEP;
  int USEnSlope;
  int USEkSlope;
  int K_BUFF;
  int N_BUFF;
  Real TEPS;
  Real EPS_L;
  Real EPS_K;
  Real EPS_N;
  Real MaxEIlmkn;
  Real MaxEIlmk;
  Real MaxEIlm;
  Real MaxEIl;
  Real MaxEHlmkn;
  Real MaxEHlmk;
  Real MaxEHlm;
  Real MaxEHl;
  Real MaxLIlmkn;
  Real MaxLIlmk;
  Real MaxLIlm;
  Real MaxLIl;
  Real MaxLHlmkn;
  Real MaxLHlmk;
  Real MaxLHlm;
  Real MaxLHl;
  Real TempMaxEIlmkn;
  int kMax;
  int nMax;

  Real *Xl;
  Real *Xlm;
  Real *Xlmk;
  Real *Xlmkn;

  ofstream outfile;

  void GetXl(int l, int ContinueSumH, int ContinueSumI);
  void GetXlm(int l, int m, int ContinueSumH, int ContinueSumI);
  void GetXlmk(int l, int m, int k, int ContinueSumH, int ContinueSumI);
  void GetXlmkn(int l, int m, int k, int n, int ContinueSumH, int ContinueSumI);


};
#endif
