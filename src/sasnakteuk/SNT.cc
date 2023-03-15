//---------------------------------------------------------------------------
//
// $Id: SNT.cc,v 1.8 2004/10/16 21:38:25 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <math.h>
#include "Globals.h"
#include "SNT.h"
#include "NRSNT.h"
#include "NRUtil.h"

SNT::SNT(const int emm, const Real rad, const Real spin,
	 const Real Omega, const Real Lambda, const Real epsilon) :
  m(emm), r(rad), a(spin), omega(Omega), lambda(Lambda), EPSILON(epsilon)
{
  EPSILON_Ext = EPSILON;

  // initialize field storage memory flag
  StoredFields = 0;

  // this is tells CalcXInf not to try to use the apoapse
  // when determining integration boundaries
  apoapse = 0.0;

  p = -666.; e = -666.;
  rp = Kerr::rplus(a);
  rm = Kerr::rminus(a);
  const Real drstardr = drsdr(r);
  const Real rstar = Kerr::rstar(r, a);
  const Real Delta = Kerr::Delta(r, a);
  const Real dr_Delta = Kerr::dr_Delta(r);

  pmomega = omega - m*a/(2.*rp);

  const Real awmm = a*omega - (Real)m;

  c0 = -12.*II*omega + lambda*(lambda + 2.) - 12.*a*omega*awmm;
  c1 =   8.*II*a*(3.*a*omega - lambda*awmm);
  c2 = -24.*II*a*awmm + 12.*a*a*(1. - 2.*awmm*awmm);
  c3 =  24.*a*a*(II*a*awmm - 1.);
  c4 = Complex(12.*a*a*a*a, 0.);

  dlmw = sqrt(2.*rp)*(8. - 12.*II*(a*m) - 4.*a*a*m*m +
		      rp*(12.*II*(a*m) + 16.*(a*m*omega - 1.) + 24.*II*omega +
			  rp*(8. - 24.*II*omega - 16.*omega*omega)));

  const Complex RInfFactor = -c0/(4.*dlmw*omega*omega);
  //
  // Code to compute derivatives of the potentials on the horizon
  //
  const Complex eta_rp = eta(rp);
  const Complex deta_rp = dr_eta(rp);
  const Real delt = rp - rm;
  const Real am = a*m;
  //
  // F deriv is easy...
  dF_rp = delt*deta_rp/(2.*rp*eta_rp);
  //
  // U deriv is a mess!
  const Complex denom = 8.*rp*rp*rp*(II*rp - am + (2.*rp*omega - II))*eta_rp;
  //
  const Complex en0 = 4.*(2.*am + II*rm)*rm*rm;
  const Complex en1 = -4*am*am*am + 4.*am*am*II*(2.*rm - 5.) +
    am*(48. - 5.*rm*rm + 2.*rm*(lambda - 16.)) -
    II*(-32. + 15.*rm*rm + 2.*(rm*(6. - lambda) - 8.*II*rm*rm*omega));
  const Complex en2 = 4.*am*am*(3.*II + 4.*omega) +
    am*(18.*rm + 56.*II*omega - 2.*(28. + lambda + 12.*II*rm*omega)) +
    II*(5.*rm*rm + 80.*II*omega + (2.*rm*(lambda - 15.) - 14.*II*rm*rm*omega)
	- 2.*(34. + lambda + rm*(24.*II*omega - 2.*II*lambda*omega)));
  const Complex en3 = am*(19. - 32.*II*omega - 15*omega*omega) -
    II*(10.*rm + 32.*omega*omega - (45. + 2.*lambda + 28.*II*rm*omega) +
	2.*omega*(48.*II + 2.*II*lambda - 8.*rm*omega));
  const Complex en4 = II*(16.*omega*omega + 34.*II*omega - 11.);

  const Complex epn1 = 2.*II*(1. - II*am)*rm*(rm - 2. - II*am);
  const Complex epn2 = 2.*II*((1. - II*am)*rm*(1. + 2.*II*omega) +
			      (rm - 2. - II*am)*(II*am - 1. + rm*
						  (2.*II*omega - 1.)));
  const Complex epn3 = 2.*II*((rm - 2. - II*am)*(1. - 2.*II*omega) +
			      (1. + 2.*II*omega)*(II*am - 1. + rm*
						  (2.*II*omega - 1.)));
  const Complex epn4 = 2.*II*(1. - 2.*II*omega)*(1. + 2.*II*omega);

  dU_rp = (eta_rp*(en0 + rp*(en1 + rp*(en2 + rp*(en3 + rp*en4)))) +
	   deta_rp*rp*(epn1 + rp*(epn2 + rp*(epn3 + rp*epn4))))/denom;
  //
  // yh[1] = r_BL
  // yh[2] = xi = log(Abs(XH))
  // yh[3] = drs_xi
  // yh[4] = phi = phase(XH)
  // yh[5] = drs_phi
  //
  Real *yh;    yh = Realvector(1, 5);
  Real *yherr; yherr = Realvector(1, 5);
  CalcXH(yh, yherr, r);
  xiH      = yh[2]; xiH_err      = yherr[2];
  drs_xiH  = yh[3]; drs_xiH_err  = yherr[3];
  psiH     = yh[4]; psiH_err     = yherr[4];
  drs_psiH = yh[5]; drs_psiH_err = yherr[5];
  free_Realvector(yh, 1, 5);
  free_Realvector(yherr, 1, 5);
  //
  XH = exp(xiH + II*psiH);
  //
  drs_XH = (drs_xiH + II*drs_psiH)*XH;
  Complex *tmp;     tmp     = Complexvector(1, 3);
  Complex *drs_tmp; drs_tmp = Complexvector(1, 3);
  tmp[1] = r;
  tmp[2] = XH;
  tmp[3] = drs_XH;
  SNT_derivs(rstar, tmp, drs_tmp);
  ddrs_XH = drs_tmp[3];
  free_Complexvector(tmp, 1, 3);
  free_Complexvector(drs_tmp, 1, 3);
  //
  dr_XH = drstardr*drs_XH;
  ddr_XH = drstardr*drstardr*ddrs_XH +
    (1./Delta)*(2.*r - dr_Delta*drstardr)*drs_XH;
  //
  Ain = CalcAin(r);
  //
  Complex *xinf;     xinf     = Complexvector(1, 3);
  Complex *xinferr;  xinferr  = Complexvector(1, 3);
  Complex *drs_xinf; drs_xinf = Complexvector(1, 3);
  CalcXInf(xinf, xinferr, r);
  XInf     = xinf[2]; XInferr     = xinferr[2];
  drs_XInf = xinf[3]; drs_XInferr = xinferr[3];
  xinf[1] = r;
  SNT_derivs(rstar, xinf, drs_xinf);
  ddrs_XInf = drs_xinf[3];
  free_Complexvector(xinf, 1, 3);
  free_Complexvector(xinferr, 1, 3);
  free_Complexvector(drs_xinf, 1, 3);
  //
  dr_XInf = drstardr*drs_XInf;
  ddr_XInf = drstardr*drstardr*ddrs_XInf +
    (1./Delta)*(2.*r - dr_Delta*drstardr)*drs_XInf;
  //
  Bin = -Ain/(4.*omega*omega);
  //
  TeukRInf = CalcR(XInf, dr_XInf, rad)*RInfFactor;
  TeukRH = CalcR(XH, dr_XH, rad);
  //
  dr_TeukRInf = Calcdr_R(XInf, dr_XInf, ddr_XInf, rad)*RInfFactor;
  dr_TeukRH = Calcdr_R(XH, dr_XH, ddr_XH, rad);
  //
  ddr_TeukRInf = Calcddr_R(TeukRInf, dr_TeukRInf, rad);
  ddr_TeukRH = Calcddr_R(TeukRH, dr_TeukRH, rad);
}

//
// This version of the constructor calculates XH quantities at the
// periapse and XInf quantites at the apoapse.
//
SNT::SNT(const int emm, const Real pee, const Real ecc, const Real spin,
	 const Real Omega, const Real Lambda, const Real epsilon) :
  m(emm), p(pee), e(ecc), a(spin), omega(Omega), lambda(Lambda),
  EPSILON(epsilon)
{
  EPSILON_Ext = EPSILON;
  
  // initialize field storage memory flag
  StoredFields = 0;

  rp = Kerr::rplus(a);
  rm = Kerr::rminus(a);
  //
  periapse = p/(1. + e);
  const Real drstardr_peri = drsdr(periapse);
  const Real rstar_peri = Kerr::rstar(periapse, a);
  const Real Delta_peri = Kerr::Delta(periapse, a);
  const Real dr_Delta_peri = Kerr::dr_Delta(periapse);
  //
  apoapse = p/(1. - e);
  const Real drstardr_apo = drsdr(apoapse);
  const Real rstar_apo = Kerr::rstar(apoapse, a);
  const Real Delta_apo = Kerr::Delta(apoapse, a);
  const Real dr_Delta_apo = Kerr::dr_Delta(apoapse);

  pmomega = omega - m*a/(2.*rp);

  const Real awmm = a*omega - (Real)m;

  c0 = -12.*II*omega + lambda*(lambda + 2.) - 12.*a*omega*awmm;
  c1 =   8.*II*a*(3.*a*omega - lambda*awmm);
  c2 = -24.*II*a*awmm + 12.*a*a*(1. - 2.*awmm*awmm);
  c3 =  24.*a*a*(II*a*awmm - 1.);
  c4 = Complex(12.*a*a*a*a, 0.);

  dlmw = sqrt(2.*rp)*(8. - 12.*II*(a*m) - 4.*a*a*m*m +
		      rp*(12.*II*(a*m) + 16.*(a*m*omega - 1.) + 24.*II*omega +
			  rp*(8. - 24.*II*omega - 16.*omega*omega)));
  //
  // Code to compute derivatives of the potentials on the horizon
  //
  const Complex eta_rp = eta(rp);
  const Complex deta_rp = dr_eta(rp);
  const Real delt = rp - rm;
  const Real am = a*m;
  //
  // F deriv is easy...
  dF_rp = delt*deta_rp/(2.*rp*eta_rp);
  //
  // U deriv is a mess!
  const Complex denom = 8.*rp*rp*rp*(II*rp - am + (2.*rp*omega - II))*eta_rp;
  //
  const Complex en0 = 4.*(2.*am + II*rm)*rm*rm;
  const Complex en1 = -4*am*am*am + 4.*am*am*II*(2.*rm - 5.) +
    am*(48. - 5.*rm*rm + 2.*rm*(lambda - 16.)) -
    II*(-32. + 15.*rm*rm + 2.*(rm*(6. - lambda) - 8.*II*rm*rm*omega));
  const Complex en2 = 4.*am*am*(3.*II + 4.*omega) +
    am*(18.*rm + 56.*II*omega - 2.*(28. + lambda + 12.*II*rm*omega)) +
    II*(5.*rm*rm + 80.*II*omega + (2.*rm*(lambda - 15.) - 14.*II*rm*rm*omega)
	- 2.*(34. + lambda + rm*(24.*II*omega - 2.*II*lambda*omega)));
  const Complex en3 = am*(19. - 32.*II*omega - 15*omega*omega) -
    II*(10.*rm + 32.*omega*omega - (45. + 2.*lambda + 28.*II*rm*omega) +
	2.*omega*(48.*II + 2.*II*lambda - 8.*rm*omega));
  const Complex en4 = II*(16.*omega*omega + 34.*II*omega - 11.);

  const Complex epn1 = 2.*II*(1. - II*am)*rm*(rm - 2. - II*am);
  const Complex epn2 = 2.*II*((1. - II*am)*rm*(1. + 2.*II*omega) +
			      (rm - 2. - II*am)*(II*am - 1. + rm*
						  (2.*II*omega - 1.)));
  const Complex epn3 = 2.*II*((rm - 2. - II*am)*(1. - 2.*II*omega) +
			      (1. + 2.*II*omega)*(II*am - 1. + rm*
						  (2.*II*omega - 1.)));
  const Complex epn4 = 2.*II*(1. - 2.*II*omega)*(1. + 2.*II*omega);

  dU_rp = (eta_rp*(en0 + rp*(en1 + rp*(en2 + rp*(en3 + rp*en4)))) +
	   deta_rp*rp*(epn1 + rp*(epn2 + rp*(epn3 + rp*epn4))))/denom;
  //
  // yh[1] = r_BL
  // yh[2] = xi = log(Abs(XH))
  // yh[3] = drs_xi
  // yh[4] = phi = phase(XH)
  // yh[5] = drs_phi
  //
  Real *yh;    yh = Realvector(1, 5);
  Real *yherr; yherr = Realvector(1, 5);
  CalcXH(yh, yherr, periapse);
  xiH      = yh[2]; xiH_err      = yherr[2];
  drs_xiH  = yh[3]; drs_xiH_err  = yherr[3];
  psiH     = yh[4]; psiH_err     = yherr[4];
  drs_psiH = yh[5]; drs_psiH_err = yherr[5];
  free_Realvector(yh, 1, 5);
  free_Realvector(yherr, 1, 5);
  //
  XH_peri = exp(xiH + II*psiH);
  XH = XH_peri;
  //
  drs_XH_peri = (drs_xiH + II*drs_psiH)*XH_peri;
  drs_XH = drs_XH_peri;
  Complex *tmp;     tmp     = Complexvector(1, 3);
  Complex *drs_tmp; drs_tmp = Complexvector(1, 3);
  tmp[1] = periapse;
  tmp[2] = XH_peri;
  tmp[3] = drs_XH_peri;
  SNT_derivs(rstar_peri, tmp, drs_tmp);
  ddrs_XH_peri = drs_tmp[3];
  free_Complexvector(tmp, 1, 3);
  free_Complexvector(drs_tmp, 1, 3);
  //
  dr_XH_peri = drstardr_peri*drs_XH_peri;
  ddr_XH_peri = drstardr_peri*drstardr_peri*ddrs_XH_peri +
    (1./Delta_peri)*(2.*periapse - dr_Delta_peri*drstardr_peri)*drs_XH_peri;
  //
  Ain = CalcAin(periapse);
  //
  Complex *xinf;     xinf     = Complexvector(1, 3);
  Complex *xinferr;  xinferr  = Complexvector(1, 3);
  Complex *drs_xinf; drs_xinf = Complexvector(1, 3);
  CalcXInf(xinf, xinferr, apoapse);
  XInf_apo     = xinf[2]; XInferr_apo     = xinferr[2];
  drs_XInf_apo = xinf[3]; drs_XInferr_apo = xinferr[3];
  xinf[1] = apoapse;
  SNT_derivs(rstar_apo, xinf, drs_xinf);
  ddrs_XInf_apo = drs_xinf[3];
  free_Complexvector(xinf, 1, 3);
  free_Complexvector(xinferr, 1, 3);
  free_Complexvector(drs_xinf, 1, 3);
  //
  dr_XInf_apo = drstardr_apo*drs_XInf_apo;
  ddr_XInf_apo = drstardr_apo*drstardr_apo*ddrs_XInf_apo +
    (1./Delta_apo)*(2.*apoapse - dr_Delta_apo*drstardr_apo)*drs_XInf_apo;
  //
  Bin = -Ain/(4.*omega*omega);
}

SNT::~SNT()
{
  if(StoredFields){
    free_Realvector(rs_stored, 1, N_stored);
    //
    free_Realvector(xiH_stored, 1, N_stored);
    free_Realvector(dxiH_stored, 1, N_stored);
    free_Realvector(ddxiH_stored, 1, N_stored);
    free_Realvector(dddxiH_stored, 1, N_stored);
    //
    free_Realvector(psiH_stored, 1, N_stored);
    free_Realvector(dpsiH_stored, 1, N_stored);
    free_Realvector(ddpsiH_stored, 1, N_stored);
    free_Realvector(dddpsiH_stored, 1, N_stored);
    //
    free_Realvector(xiInf_stored, 1, N_stored);
    free_Realvector(dxiInf_stored, 1, N_stored);
    free_Realvector(ddxiInf_stored, 1, N_stored);
    free_Realvector(dddxiInf_stored, 1, N_stored);
    //
    free_Realvector(psiInf_stored, 1, N_stored);
    free_Realvector(dpsiInf_stored, 1, N_stored);
    free_Realvector(ddpsiInf_stored, 1, N_stored);
    free_Realvector(dddpsiInf_stored, 1, N_stored);
  }
}
