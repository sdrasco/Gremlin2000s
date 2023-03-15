//---------------------------------------------------------------------------
//
// $Id: LB.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <math.h>
#include "Globals.h"
#include "LB.h"
#include "NRLB.h"
#include "NRUtil.h"

LB::LB(const Real arr, const Real ayy) : r(arr), a(ayy)
{
  gkg = new GKG;

  const Real isco_ret = Kerr::isco_ret(a);

  if (r < isco_ret) {
    const Real tolx = 1.e-14;
    const Real tolf = 1.e-14;
    //
    Real *x; x = Realvector(1, 3);
    x[1] = Kerr::Eeqpro(r, a);
    x[2] = Kerr::Lzeqpro(r, a);
    x[3] = 0.;
    //
    mnewt(200, x, 3, tolx, tolf);
    //
    E = x[1];
    Lz = x[2];
    Q = x[3];
    cosiota = Lz/sqrt(Lz*Lz + Q);
    free_Realvector(x, 1, 3);
  } else {
    E = Kerr::Eeqret(r, a);
    Lz = Kerr::Lzeqret(r, a);
    Q = 0.;
    cosiota = -1.;
  }
}

LB::~LB()
{
  delete gkg;
}

void LB::mnewt(int ntrial, Real x[], int n, Real tolx, Real tolf)
{
  int k, i, *indx;
  Real errx, errf, d, *fvec, **fjac, *p;
  
  indx = ivector(1, n);
  p = Realvector(1, n);
  fvec = Realvector(1, n);
  fjac = Realmatrix(1, n, 1, n);
  for (k = 1; k <= ntrial; k++) {
    func_and_jac(x, n, fvec, fjac);
    errf = 0.0;
    for (i = 1; i <= n; i++) errf += fabs(fvec[i]);
    if (errf <= tolf) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
    for (i = 1; i <= n; i++) p[i] = -fvec[i];
    ludcmp(fjac, n, indx, &d);
    lubksb(fjac, n, indx, p);
    errx = 0.0;
    for (i = 1; i <= n; i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
    }
    if (errx <= tolx) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
  }
  free_Realmatrix(fjac, 1, n, 1, n);
  free_Realvector(fvec, 1, n);
  free_Realvector(p, 1, n);
  free_ivector(indx, 1, n);
  return;
}

void LB::func_and_jac(Real *x, int n, Real *fvec, Real **fjac)
{
  if (n == 3) {
    const Real E = x[1];
    const Real Lz = x[2];
    const Real Q = x[3];
    
    fvec[1] = gkg->RFunc(r, a, 0., E, Lz, Q);
    fvec[2] = gkg->dr_RFunc(r, a, 0., E, Lz, Q);
    fvec[3] = gkg->ddr_RFunc(r, a, 0., E, Lz, Q);
    //
    // No separate calls for the Jacobian elements ...
    // they're fairly simple, though.
    //
    const Real tmp1 = r*r + a*a;
    const Real tmp2 = Lz - a*E;
    
    fjac[1][1] = 2.*r*(E*r*tmp1 - 2.*a*tmp2);
    fjac[1][2] = 2.*r*(2.*tmp2 - Lz*r);
    fjac[1][3] = 2.*r  - tmp1;
    fjac[2][1] = 4.*(E*r*(a*a + 2.*r*r) - a*tmp2);
    fjac[2][2] = 4.*(tmp2 - Lz*r);
    fjac[2][3] = 2.*(1. - r);
    fjac[3][1] = 4.*E*(a*a + 6.*r*r);
    fjac[3][2] = -4.*Lz;
    fjac[3][3] = -2.;
  } else if (n == 2) {
    const Real E = 1.;
    const Real Lz = x[1];
    const Real Q = x[2];
    
    fvec[1] = gkg->RFunc(r, a, 0., E, Lz, Q);
    fvec[2] = gkg->dr_RFunc(r, a, 0., E, Lz, Q);
    //
    // No separate calls for the Jacobian elements ...
    // they're fairly simple, though.
    //
    fjac[1][1] = -2.*(2.*a + Lz*(r - 2.))*r;
    fjac[1][2] = -Kerr::Delta(r, a);
    fjac[2][1] = -4.*(a + Lz*(r - 1.));
    fjac[2][2] = 2.*(1. - r);
  } else {
    cerr << "Case n = " << n << " in LBFJ not understood ...";
    exit(0);
  }
}

