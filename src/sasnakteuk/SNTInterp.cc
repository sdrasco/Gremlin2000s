//---------------------------------------------------------------------------
//
// $Id: SNTInterp.cc,v 1.9 2004/06/16 18:59:09 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// This code uses spline interpolation to compute the fields R^H 
// and R^Inf (and their first two r-derivatives) for arbitrary r.
// 
// To use it, you must first call PrepSpline(N) which stores the results of 
// calling CalcFields at N evenly spaced intervals which enclose 
// orbital range in r.  After that, you can call the interpolation
// routine SplineFields(r,DoH) in exactly the same way that you would 
// call CalcFields(r,DoH) 
//
// Note: DoH = 1 means compute the H-fields
//       otherwise we compute the Z-fields 
//  
// Steve Drasco, 29 May 2004
//
// Modified somewhat heavily by Scott Hughes, 30 May 2004:
//
//     1. Rather than interpolating in R(r), it interpolates in X(r*).
//         a. More accurately, it interpolates in xi and psi, defined
//            via X = exp(xi + I psi).  In other words, xi is the log
//            amplitude, psi is the phase.  These variables are very
//            smooth.
//         b. It uses rstar as the radial coordinate for the spline
//            rather than r.  Things are very smooth in rstar.
//     2. Second derivatives are computed numerically with spline()
//        rather than analytically.  I tried both, and this seems to
//        be more accurate.  Go figure.
//     3. Changed domain over which spline is constructed to exactly
//        match the radial motion, rather 10% of padding on each side.
//        The padding seemed to be causing errors, particularly near
//        periapse.
//
#include <math.h>
#include "Globals.h"
#include "SNT.h"

void SNT::PrepSpline(int N) {
  // make sure we have the right kind of object
  if (e == -666.)
    Die("The wrong constructor was called ... waaaahhh");

  // set the number of points
  N_stored = N;

  // allocate memory
  rs_stored = Realvector(1, N);
  //
  xiH_stored = Realvector(1, N);
  dxiH_stored = Realvector(1, N);
  ddxiH_stored = Realvector(1, N);
  dddxiH_stored = Realvector(1, N);
  //
  psiH_stored = Realvector(1, N);
  dpsiH_stored = Realvector(1, N);
  ddpsiH_stored = Realvector(1, N);
  dddpsiH_stored = Realvector(1, N);
  //
  xiInf_stored = Realvector(1, N);
  dxiInf_stored = Realvector(1, N);
  ddxiInf_stored = Realvector(1, N);
  dddxiInf_stored = Realvector(1, N);
  //
  psiInf_stored = Realvector(1, N);
  dpsiInf_stored = Realvector(1, N);
  ddpsiInf_stored = Realvector(1, N);
  dddpsiInf_stored = Realvector(1, N);
  //
  // calculate the fields at N points
  //
  Real rs_start = Kerr::rstar(p/(1.0 + e), a);
  //
  Real rs_finish = Kerr::rstar(p/(1.0 - e), a);
  Real drs = (rs_finish - rs_start)/Real(N - 1);
  Real rad;
  //
  Complex dlnXH, dlnXInf;
  for(int i = 1; i <= N; i++){
    rs_stored[i] = rs_start + Real(i - 1)*drs;
    //
    // Quick bisection to get r at this rstar.
    //
    Real rlo = 0.999*p/(1.0 + e), rhi = 1.001*p/(1.0 - e);
    while ((rhi - rlo)/rlo > 4.e-16) {
      if (Kerr::rstar(0.5*(rhi + rlo), a) > rs_stored[i]) {
	rhi = 0.5*(rhi + rlo);
      } else {
	rlo = 0.5*(rhi + rlo);
      }
    }
    rad = 0.5*(rhi + rlo);
    //
    CalcFields(rad, 1);
    xiH_stored[i] = log(abs(XH));
    if (i > 1)
      psiH_stored[i] = my_atan2(XH.imag(), XH.real(), psiH_stored[i - 1]);
    else
      psiH_stored[i] = my_atan2(XH.imag(), XH.real(), 0.);
    dlnXH = drs_XH/XH;
    dxiH_stored[i] = dlnXH.real();
    dpsiH_stored[i] = dlnXH.imag();
    //
    CalcFields(rad, 0);
    xiInf_stored[i] = log(abs(XInf));
    if (i > 1)
      psiInf_stored[i] = my_atan2(XInf.imag(), XInf.real(),
				  psiInf_stored[i - 1]);
    else
      psiInf_stored[i] = my_atan2(XInf.imag(), XInf.real(), 0.);
    dlnXInf = drs_XInf/XInf;
    dxiInf_stored[i] = dlnXInf.real();
    dpsiInf_stored[i] = dlnXInf.imag();
  }
  spline(rs_stored, xiH_stored, N, 1.e30, 1.e30, ddxiH_stored);
  spline(rs_stored, psiH_stored, N, 1.e30, 1.e30, ddpsiH_stored);
  spline(rs_stored, dxiH_stored, N, 1.e30, 1.e30, dddxiH_stored);
  spline(rs_stored, dpsiH_stored, N, 1.e30, 1.e30, dddpsiH_stored);
  //
  spline(rs_stored, xiInf_stored, N, 1.e30, 1.e30, ddxiInf_stored);
  spline(rs_stored, psiInf_stored, N, 1.e30, 1.e30, ddpsiInf_stored);
  spline(rs_stored, dxiInf_stored, N, 1.e30, 1.e30, dddxiInf_stored);
  spline(rs_stored, dpsiInf_stored, N, 1.e30, 1.e30, dddpsiInf_stored);
  //
  // set flag for freeing memory in the destructor
  StoredFields = 1;
}

Real SNT::my_atan2(const Real y, const Real x, const Real prev)
{
  Real ans = atan2(y, x);
  Real factor = 0.95*(2.*M_PI);

  if (fabs(ans - prev) > factor) {
    //
    // Need to do some smoothing!
    if (ans - prev < 0) {
      //
      // ans is too small: add 2Pi repeatedly.
      while (ans - prev < -factor)
	ans += 2.*M_PI;
    } else {
      //
      // ans is too big: subtract 2Pi repeatedly.
      while (ans - prev > factor)
	ans -= 2.*M_PI;
    }
  }
  return(ans);
}

void SNT::SplineFields(Real rad, int DoH){
  //
  // make sure we have the right kind of object
  if (e == -666.)
    Die("The wrong constructor was called ... waaaahhh");
  //
  // make sure we stored the fields 
  if (StoredFields == 0) 
    Die("You must use SNT:PrepSpline before you use SNT::SplineFields");
  //
  Real radstar = Kerr::rstar(rad, a);
  Real Delta = Kerr::Delta(rad, a);
  Real dr_Delta = Kerr::dr_Delta(rad);
  Real xiH, psiH, xiInf, psiInf;
  Real dxiH, dpsiH, dxiInf, dpsiInf;
  Real drstardr = drsdr(rad);
  Complex *X, *dX;
  X = Complexvector(1, 3); dX = Complexvector(1, 3);
  //
  // compute the fields
  if (DoH) {
    splint(rs_stored, xiH_stored, ddxiH_stored, N_stored, radstar, &xiH);
    splint(rs_stored, psiH_stored, ddpsiH_stored, N_stored, radstar, &psiH);
    splint(rs_stored, dxiH_stored, dddxiH_stored, N_stored, radstar, &dxiH);
    splint(rs_stored, dpsiH_stored, dddpsiH_stored, N_stored, radstar, &dpsiH);
    XH = exp(xiH + II*psiH);
    drs_XH = (dxiH + II*dpsiH)*XH;
    dr_XH = drs_XH*drstardr;
    //
    // Need second derivatives to get Teukolsky functions!
    X[1] = rad; X[2] = XH; X[3] = drs_XH;
    SNT_derivs(radstar, X, dX);
    ddrs_XH = dX[3];
    ddr_XH = drstardr*drstardr*ddrs_XH +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XH;
    TeukRH = CalcR(XH, dr_XH, rad);
    dr_TeukRH = Calcdr_R(XH, dr_XH, ddr_XH, rad);
    ddr_TeukRH = Calcddr_R(TeukRH, dr_TeukRH, rad);
  } else {
    //
    const Complex RInfFactor = -c0/(4.*dlmw*omega*omega);
    //
    splint(rs_stored, xiInf_stored, ddxiInf_stored, N_stored, radstar, &xiInf);
    splint(rs_stored, psiInf_stored, ddpsiInf_stored, N_stored, radstar,
	   &psiInf);
    splint(rs_stored, dxiInf_stored, dddxiInf_stored, N_stored, radstar,
	   &dxiInf);
    splint(rs_stored, dpsiInf_stored, dddpsiInf_stored, N_stored, radstar,
	   &dpsiInf);
    XInf = exp(xiInf + II*psiInf);
    drs_XInf = (dxiInf + II*dpsiInf)*XInf;
    dr_XInf = drs_XInf*drstardr;
    //
    // Need second derivatives to get Teukolsky functions!
    X[1] = rad; X[2] = XInf; X[3] = drs_XInf;
    SNT_derivs(radstar, X, dX);
    ddrs_XInf = dX[3];
    ddr_XInf = drstardr*drstardr*ddrs_XInf +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XInf;
    TeukRInf = CalcR(XInf, dr_XInf, rad)*RInfFactor;
    dr_TeukRInf = Calcdr_R(XInf, dr_XInf, ddr_XInf, rad)*RInfFactor;
    ddr_TeukRInf = Calcddr_R(TeukRInf, dr_TeukRInf, rad);
  }
  free_Complexvector(X, 1, 3); free_Complexvector(dX, 1, 3);
}
