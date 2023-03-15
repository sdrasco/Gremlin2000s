//---------------------------------------------------------------------------
//
// $Id: SNTSasNak.cc,v 1.12 2004/11/08 20:29:31 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Code for the Sasaki-Nakamura equation derivatives and solver.
//
// Scott Hughes, 8 January 1999
//
#include <math.h>
#include <iostream.h>
#include "Globals.h"
#include "SNT.h"
#include "NRUtil.h"
#include "NRSNT.h"

#define ROUT_OMEGA 35.

void SNT::CalcXInf(Complex XInf[], Complex XInferr[], const Real rad)
{
  int i = 1, itmax = 15;
  Real rout, rsout;
  Real radstar = Kerr::rstar(rad, a);
  Complex **XInfest;    XInfest = Complexmatrix(1, itmax, 1, 3);
  Complex **d;                d = Complexmatrix(1, itmax, 1, itmax);
  Real *x;                    x = Realvector(1, itmax);

// start of steve's changes
  rout = Max(ROUT_OMEGA/fabs(omega), 2.0*rad); // original code
  if(rout < apoapse) rout = 2*apoapse; // new
// end of steve's changes
  do {
    rsout = Kerr::rstar(rout, a);
    XInfInitCond(XInfest[i], rout, rsout);
    IntegrateODE(XInfest[i], rsout, radstar, EPSILON_Int);
    rzextr(i, pow(0.5, i-1), XInfest[i], XInf, XInferr, d, x, 3);
    i++;
    rout *= 2.0;
  } while((fabs(XInferr[2].real()/XInf[2].real()) > EPSILON_Ext ||
	   fabs(XInferr[2].imag()/XInf[2].imag()) > EPSILON_Ext ||
	   fabs(XInferr[3].real()/XInf[3].real()) > EPSILON_Ext ||
	   fabs(XInferr[3].imag()/XInf[3].imag()) > EPSILON_Ext) &&
	   i < itmax);

  free_Realvector(x, 1, itmax);
  free_Complexmatrix(d, 1, itmax, 1, itmax);
  free_Complexmatrix(XInfest, 1, itmax, 1, 3);
}

//
// This routine assumes that CalcXH has been run; it picks up
// from where it left off.
//
Complex SNT::CalcAin(const Real rad)
{
  int i = 1, itmax = 15, j;
  Real rin = rad, rout, rout_orig;
  Complex *X;  X = Complexvector(1, 3);
  Complex **Ainest;  Ainest = Complexmatrix(1, itmax, 1, 1);
  Complex *Ainextrp;  Ainextrp = Complexvector(1, 1);
  Complex *dAinextrp;  dAinextrp = Complexvector(1, 1);
  //
  // Start at the field point, using XH as the initial condition.
  //
  X[1] = rin;
  X[2] = XH;
  X[3] = drs_XH;
  //
  // Used to pull Ain from integrated SN function
  //
  Real rsin, rsout, drdrs;
  Complex PP, PPb, dr_PP, dr_PPb, eiwrs;
  Complex numer, denom;
  //
  // Used for the extrapolation
  //
  Complex **d;  d = Complexmatrix(1, itmax, 1, itmax);
  Real *x;  x = Realvector(1, itmax);
  Real scale;
// start of steve's changes
  rout_orig = Max(ROUT_OMEGA/fabs(omega), 2.0*rad); // original code
  if(rout_orig < apoapse) rout_orig = 2*apoapse; // new
// end of steve's changes
  rsin = Kerr::rstar(rin, a);
  rsout = Kerr::rstar(rout_orig, a);
  do {
    IntegrateODE(X, rsin, rsout, EPSILON_Int);
    rout = X[1].real();
    eiwrs = exp(II*omega*rsout);

    PP = P(rout);
    PPb = conj(PP);
    dr_PP = dr_P(rout);
    dr_PPb = conj(dr_PP);
    drdrs = 1./drsdr(rout);

    numer = -X[3]*PPb + X[2]*(II*omega*PPb + drdrs*dr_PPb);
    denom =  2.*II*omega*PP*PPb + drdrs*(PP*dr_PPb - PPb*dr_PP);
    Ainest[i][1] = eiwrs*numer/denom;

    rzextr(i, pow(rout_orig/rout, i-1), Ainest[i], Ainextrp,
	   dAinextrp, d, x, 1);

    rsin = rsout;
    rsout *= 2.0;
    i++;
  } while((fabs(dAinextrp[1].real()/Ainextrp[1].real()) > EPSILON_Ext ||
	   fabs(dAinextrp[1].imag()/Ainextrp[1].imag()) > EPSILON_Ext) &&
	  i < itmax);
  free_Complexvector(X, 1, 3);
  free_Complexmatrix(Ainest, 1, itmax, 1, 1);
  free_Complexmatrix(d, 1, itmax, 1, itmax);
  free_Realvector(x, 1, itmax);

  const Complex ans = Ainextrp[1];

  free_Complexvector(Ainextrp, 1, 1);
  free_Complexvector(dAinextrp, 1, 1);

  return(ans);
}

void SNT::CalcXH(Real YH[], Real YHerr[], const Real rad)
{
  int i, itmax = 9;
  Real rsin, rsout, eps;
  Real rin, rout;
  Real drstar_dr;
  Real **YHest; YHest = Realmatrix(1, itmax, 1, 5);
  //
  // Used for the polynomial extrapolation
  //
  Real **d; d = Realmatrix(1, itmax, 1, itmax);
  Real *x;  x = Realvector(1, itmax);
  //
  const Real scale = 0.002*(rp/rad);
  rsout = Kerr::rstar(rad, a);
  i = 1;
  do {
    eps = EPSILON/pow(2, (Real)(i-2));
    rin = rp + eps;
    rsin = Kerr::rstar(rin, a);
    XHInitCond(YHest[i], eps, rsin);
    IntegrateODE(YHest[i], rsin, rsout, 10.*EPSILON_Int);
    rzextr(i, pow(scale, (Real)(i-1)), YHest[i], YH, YHerr, d, x, 5);
    i++;
  } while((fabs(YHerr[1]/YH[1]) > EPSILON_Ext ||
	   fabs(YHerr[2]/YH[2]) > EPSILON_Ext ||
	   fabs(YHerr[3]/YH[3]) > EPSILON_Ext ||
	   fabs(YHerr[4]/YH[4]) > EPSILON_Ext ||
	   fabs(YHerr[5]/YH[5]) > EPSILON_Ext) &&
	   i < itmax);
  free_Realvector(x, 1, itmax);
  free_Realmatrix(d, 1, itmax, 1, itmax);
  free_Realmatrix(YHest, 1, itmax, 1, 5);
}

void SNT::CalcFields(const Real rad, int DoH)
{
  if (e == -666.) Die("The wrong constructor was called ... waaaahhh");
  const Real Delta = Kerr::Delta(rad, a);
  const Real dr_Delta = Kerr::dr_Delta(rad);
  const Real drstardr = drsdr(rad);
  //
  Complex *X;  X = Complexvector(1, 3);
  Complex *tmp, *drs_tmp;
  Real rsin, rsout;
  rsout = Kerr::rstar(rad, a);

  if (DoH) {
    //
    // Calculate XH.  Start at periapse, using XH_peri as the
    // initial condition.
    //
    X[1] = periapse;
    X[2] = XH_peri;
    X[3] = drs_XH_peri;
    //
    rsin = Kerr::rstar(periapse, a);
    //
    if (rsout - rsin > 1e-13)
      IntegrateODE(X, rsin, rsout, EPSILON_Int);
    XH = X[2];
    drs_XH = X[3];
    //
    tmp     = Complexvector(1, 3);
    drs_tmp = Complexvector(1, 3);
    tmp[1] = rad;
    tmp[2] = XH;
    tmp[3] = drs_XH;
    SNT_derivs(rsout, tmp, drs_tmp);
    ddrs_XH = drs_tmp[3];
    free_Complexvector(tmp, 1, 3);
    free_Complexvector(drs_tmp, 1, 3);
    //
    dr_XH = drstardr*drs_XH;
    ddr_XH = drstardr*drstardr*ddrs_XH +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XH;
  } else {
    //
    // Calculate XInf.  Start at apoapse, using XH_apo as the
    // initial condtion.
    // 
    X[1] = apoapse;
    X[2] = XInf_apo;
    X[3] = drs_XInf_apo;
    //
    rsin = Kerr::rstar(apoapse, a);
    //
    if (rsin - rsout > 1e-13)
      IntegrateODE(X, rsin, rsout, EPSILON_Int);
    XInf = X[2];
    drs_XInf = X[3];
    //
    tmp     = Complexvector(1, 3);
    drs_tmp = Complexvector(1, 3);
    tmp[1] = rad;
    tmp[2] = XInf;
    tmp[3] = drs_XInf;
    SNT_derivs(rsout, tmp, drs_tmp);
    ddrs_XInf = drs_tmp[3];
    free_Complexvector(tmp, 1, 3);
    free_Complexvector(drs_tmp, 1, 3);
    //
    dr_XInf = drstardr*drs_XInf;
    ddr_XInf = drstardr*drstardr*ddrs_XInf +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XInf;
  }
  free_Complexvector(X, 1, 3);
}
//
// This computes the R fields as well as the X fields
//
void SNT::CalcRXFields(const Real rad, int DoH)
{
  if (e == -666.) Die("The wrong constructor was called ... waaaahhh");
  const Real Delta = Kerr::Delta(rad, a);
  const Real dr_Delta = Kerr::dr_Delta(rad);
  const Real drstardr = drsdr(rad);
  //
  Complex *X;  X = Complexvector(1, 3);
  Complex *tmp, *drs_tmp;
  Real rsin, rsout;
  rsout = Kerr::rstar(rad, a);

  if (DoH) {
    //
    // Calculate H fields.  Start at periapse, using XH_peri as the
    // initial condition.
    //
    X[1] = periapse;
    X[2] = XH_peri;
    X[3] = drs_XH_peri;
    //
    rsin = Kerr::rstar(periapse, a);
    //
    if (rsout - rsin > 1e-13)
      IntegrateODE(X, rsin, rsout, EPSILON_Int);
    XH = X[2];
    drs_XH = X[3];
    //
    tmp     = Complexvector(1, 3);
    drs_tmp = Complexvector(1, 3);
    tmp[1] = rad;
    tmp[2] = XH;
    tmp[3] = drs_XH;
    SNT_derivs(rsout, tmp, drs_tmp);
    ddrs_XH = drs_tmp[3];
    free_Complexvector(tmp, 1, 3);
    free_Complexvector(drs_tmp, 1, 3);
    //
    dr_XH = drstardr*drs_XH;
    ddr_XH = drstardr*drstardr*ddrs_XH +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XH;
    //
    TeukRH = CalcR(XH, dr_XH, rad);
    dr_TeukRH = Calcdr_R(XH, dr_XH, ddr_XH, rad);
    ddr_TeukRH = Calcddr_R(TeukRH, dr_TeukRH, rad);
  } else {
    //
    // Calculate Inf fields.  Start at apoapse, using XH_apo as the
    // initial condition.
    //
    X[1] = apoapse;
    X[2] = XInf_apo;
    X[3] = drs_XInf_apo;
    const Complex RInfFactor = -c0/(4.*dlmw*omega*omega);
    //
    rsin = Kerr::rstar(apoapse, a);
    //
    if (rsin - rsout > 1e-13)
      IntegrateODE(X, rsin, rsout, EPSILON_Int);
    XInf = X[2];
    drs_XInf = X[3];
    //
    tmp     = Complexvector(1, 3);
    drs_tmp = Complexvector(1, 3);
    tmp[1] = rad;
    tmp[2] = XInf;
    tmp[3] = drs_XInf;
    SNT_derivs(rsout, tmp, drs_tmp);
    ddrs_XInf = drs_tmp[3];
    free_Complexvector(tmp, 1, 3);
    free_Complexvector(drs_tmp, 1, 3);
    //
    dr_XInf = drstardr*drs_XInf;
    ddr_XInf = drstardr*drstardr*ddrs_XInf +
      (1./Delta)*(2.*rad - dr_Delta*drstardr)*drs_XInf;
    //
    TeukRInf = CalcR(XInf, dr_XInf, rad)*RInfFactor;
    dr_TeukRInf = Calcdr_R(XInf, dr_XInf, ddr_XInf, rad)*RInfFactor;
    ddr_TeukRInf = Calcddr_R(TeukRInf, dr_TeukRInf, rad);
  }
  free_Complexvector(X, 1, 3);
}

void SNT::IntegrateODE(Real y[], const Real rs1, const Real rs2,
		       const Real tol)
{
  const Real drs = fabs(rs1 - rs2);

  odeint(y, 5, rs1, rs2, tol, drs, tol/100.);
}

void SNT::IntegrateODE(Complex X[], const Real rs1, const Real rs2,
		       const Real tol)
{
  const Real drs = fabs(rs1 - rs2);

  odeint(X, 3, rs1, rs2, tol, drs, tol/100.);
}

void SNT::SNT_derivs(const Real rs, Real Y[], Real drs_Y[])
{
  const Real Fr = F(Y[1]).real();
  const Real Fi = F(Y[1]).imag();
  const Real Ur = U(Y[1]).real();
  const Real Ui = U(Y[1]).imag();
  const Real Delta = Kerr::Delta(Y[1], a);
  const Real r2pa2 = Y[1]*Y[1] + a*a;
  //
  drs_Y[1] = Delta/r2pa2;
  //
  // Y[2] = xi    Y[3] = dxi/drs
  // Y[4] = psi   Y[5] = dpsi/drs
  //
  drs_Y[2] = Y[3];
  drs_Y[4] = Y[5];
  //
  drs_Y[3] = Y[5]*Y[5] - Y[3]*Y[3] + Fr*Y[3] - Fi*Y[5] + Ur;
  drs_Y[5] = -2.*Y[3]*Y[5] + Fi*Y[3] + Fr*Y[5] + Ui;
}

void SNT::SNT_derivs(const Real rs, Complex X[], Complex drs_X[])
{
  const Real rad = X[1].real();
  const Real Delta = Kerr::Delta(rad, a);
  const Real r2pa2 = rad*rad + a*a;

  drs_X[1] = Delta/r2pa2;
  drs_X[2] = X[3];
  drs_X[3] = U(rad)*X[2] + F(rad)*X[3];
}

void SNT::XInfInitCond(Complex XInf_Init[], const Real rout, const Real rsout)
{
  const Complex PPb = conj(P(rout));
  const Complex eiwrs = exp(II*omega*rsout);

  XInf_Init[1] = rout;
  XInf_Init[2] = eiwrs*PPb;
  XInf_Init[3] = eiwrs*(II*omega*PPb + conj(dr_P(rout))/drsdr(rout));
}

void SNT::XHInitCond(Real YH_Init[], const Real deltar, const Real rsin)
{
  const Real d = rp - rm;
  const Complex A = 4.*rp*rp*(II*dU_rp + pmomega*dF_rp)/
    (d*(II*d + 4.*pmomega*rp));
  const Complex eiprs = exp(-II*pmomega*rsin);
  //
  // First, the complex version.
  //
  const Complex y = (1. + deltar*A)*eiprs; 
  const Complex dy = (-II*pmomega*(1. + deltar*A) + A*deltar*d/(2.*rp))*eiprs;
  const Complex dlny = dy/y;
  //
  // The B-L radius ...
  //
  YH_Init[1] = deltar + rp;
  //
  // ... log amplitude and phase ...
  //
  YH_Init[2] = log(abs(y));
  YH_Init[4] = atan2(y.imag(), y.real());
  //
  // ... and derivatives ...
  //
  YH_Init[3] = dlny.real();
  YH_Init[5] = dlny.imag();
}

Complex SNT::P(const Real rad)
{
  const Real x = 1./(omega*rad);
  const Real lp2 = lambda + 2.;
  const Real amw = a*m*omega;
  const Complex Acof = -II*(lp2 + 2.*amw)/2.;
  const Complex Bcof = (-lp2*lp2 + lp2*(2. - 4.*amw) +
 			4.*(amw + 3.*omega*II - amw*(amw + 2.*omega*II)))/8.;
  const Complex Ccof = II*(-4.*amw - Bcof*(lambda - 4. + 2.*amw +
					   8.*II*omega) +
			   omega*(-12.*omega -
				  II*lambda*(lp2 + 2.*amw) + 
				  a*a*omega*((Real)(m*m) + lambda -
					     3. + 2.*amw)))/6.;

  return(1. + x*(Acof + x*(Bcof + x*Ccof)));
}

Complex SNT::dr_P(const Real rad)
{
  const Real x = 1./(omega*rad);
  const Real lp2 = lambda + 2.;
  const Real amw = a*m*omega;
  const Complex Acof = -II*(lp2 + 2.*amw)/2.;
  const Complex Bcof = (-lp2*lp2 + lp2*(2. - 4.*amw) +
  			4.*(amw + 3.*omega*II - amw*(amw + 2.*omega*II)))/8.;
  const Complex Ccof = II*(-4.*amw - Bcof*(lambda - 4. + 2.*amw +
					   8.*II*omega) +
			   omega*(-12.*omega -
				  II*lambda*(lp2 + 2.*amw) + 
				  a*a*omega*((Real)(m*m) + lambda -
					     3. + 2.*amw)))/6.;
  
  return(-omega*x*x*(Acof + x*(2.*Bcof + 3.*x*Ccof)));
}

Real SNT::drsdr(const Real rad)
{
  return((rad*rad + a*a)/Kerr::Delta(rad, a));
}
