//---------------------------------------------------------------------------
//
// $Id: SNT.h,v 1.8 2004/06/16 18:59:09 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// SNT --- potentials and routines for solving the Sasaki-Nakamura
// equation, and then converting the Sasaki-Nakamura solutions to the
// Teukolsky solution.
//
// Scott Hughes, started 7 January 1999.
//
// Major modifications begun 31 December 2002: goal is to change the
// integrator so that it integrates the phase and amplitude of the SN
// equation separately, rather than as a mess of complex variables.
//
// Wrapping up the modifications, 21 Feb 2003: after extensive
// testing, I believe I've converged to a set of equations that
// perform rather well.  We use the phase and amplitude formulation
// for the near horizon region, and thus to calculate XH.  We then use
// the old complex number integrator for computing XInf and Ain.  Note
// that we also integrate all quantities in rstar, rather than the
// Boyer-Lindquist coordinate r.  In order to have r_BL available for
// evaluating the potentials, we introduce additional fields to
// integrate r(rstar).
//
#ifndef _SNT_H
#define _SNT_H
#include "Globals.h"
#include "NRUtil.h"

#define EPSILON_Int 1.e-14  // Epsilon used in all integration.

// this is in numrec
void splint(Real xa[], Real ya[], Real y2a[], int n, Real x, Real *y);
void spline(Real x[], Real y[], int n, Real yp1, Real ypn, Real y2[]);

class SNT {
public:
  SNT(const int emm, const Real rad,
      const Real spin, const Real Omega,
      const Real Lambda, const Real epsilon);
  SNT(const int emm, const Real pee, const Real ecc,
      const Real spin, const Real Omega,
      const Real Lambda, const Real epsilon);
  ~SNT();
  //
  // Leaving almost everything public now for debugging purposes.
  // A lot of this will be "privatized" later.
  //
  Complex Ain, Bin;
  Complex TeukRH, TeukRInf;
  Complex dr_TeukRH, dr_TeukRInf;
  Complex ddr_TeukRH, ddr_TeukRInf;
  //
  // Note there are no error variables associated with XH; we
  // keep them associated with xi and psi.
  //
  Real xiH, drs_xiH;
  Real psiH, drs_psiH;
  Real xiH_err, drs_xiH_err;
  Real psiH_err, drs_psiH_err;
  Complex XH, XInf, XInferr;
  Complex drs_XH, drs_XInf, drs_XInferr;
  Complex ddrs_XH, ddrs_XInf;
  Complex dr_XH, dr_XHerr, dr_XInf, dr_XInferr;
  Complex ddr_XH, ddr_XInf;
  //
  // For eccentric orbits, the quantities that we compute in the
  // constructor are XH (and associated quantities) at periapse,
  // and XInf (and associated quantities) at apoapse.
  //
  Complex XH_peri, XInf_apo, XInferr_apo;
  Complex drs_XH_peri, drs_XInf_apo, drs_XInferr_apo;
  Complex ddrs_XH_peri, ddrs_XInf_apo;
  Complex dr_XH_peri, dr_XHerr_peri, dr_XInf_apo, dr_XInferr_apo;
  Complex ddr_XH_peri, ddr_XInf_apo;
  //
  // This is for spline interpolation of the fields (IE orbits)
  // 
  int StoredFields, N_stored;
  Real *rs_stored;
  Real *xiH_stored, *ddxiH_stored;
  Real *psiH_stored, *ddpsiH_stored;
  Real *xiInf_stored, *ddxiInf_stored;
  Real *psiInf_stored, *ddpsiInf_stored;
  Real *dxiH_stored, *dddxiH_stored;
  Real *dpsiH_stored, *dddpsiH_stored;
  Real *dxiInf_stored, *dddxiInf_stored;
  Real *dpsiInf_stored, *dddpsiInf_stored;
  void PrepSpline(int N);
  void SplineFields(Real rad, int DoH);
  //
  // These must be left public for lots of stuff in the
  // Green's function.
  //
  inline Real K(const Real rad)
    {
      return((rad*rad + a*a)*omega - m*a);
    };

  inline Real dr_K(const Real rad)
    {
      return(2.*rad*omega);
    };

  inline Real ddr_K()
    {
      return(2.*omega);
    };
  //
  // This next routine is used to compute all relevant field quantities
  // (including RH, RInf, XH, XInf, and derivatives) at a radius between
  // periapse and apoapse.
  //
  void CalcFields(const Real rad, int DoH); // computes X fields only
  void CalcRXFields(const Real rad, int DoH); // computes R & X fields

private:
  int m;
  int imax;
  //
  // These variables are only used when we are working with
  // eccentric orbits.
  //
  Real p, e, periapse, apoapse;
  //
  Real r, a, rp, rm, rswitch, omega, lambda, pmomega;
  Complex dF_rp, dU_rp; // 1st derivs of potentials at horizon.
  Real EPSILON, EPSILON_Ext;
  Complex dlmw;

  Complex c0, c1, c2, c3, c4;

  inline Complex eta(const Real rad)
    {
      const Real x = 1./rad;

      return(c0 + x*(c1 + x*(c2 + x*(c3 + x*c4))));
    };

  inline Complex dr_eta(const Real rad)
    {
      const Real x = 1./rad;

      return(-x*x*(c1 + x*(2.*c2 + x*(3.*c3 + x*4.*c4))));
    };

  //
  // These two functions are given in the simplified form found by
  // Mathematica --- quite a bit simpler than their "canonical" form.
  //
  inline Real G(const Real rad)
    {
      const Real r2pa2 = rad*rad + a*a;
      
      const Real numer = a*a*(2. - rad) - rad*rad*rad;
      const Real denom = r2pa2*r2pa2;
      
      return(numer/denom);
    };

  inline Real dr_G(const Real rad)
    {
      const Real r2pa2 = rad*rad + a*a;
      
      const Real numer = rad*rad*rad*rad - a*a*(a*a + 8.*rad);
      const Real denom = r2pa2*r2pa2*r2pa2;
      
      return(numer/denom);
    };

  inline Complex V(const Real rad)
    {
      const Real Kay = K(rad);
      
      const Complex V_1 = -Kay*(Kay + 4.*II*(rad - 1.))/
	Kerr::Delta(rad, a);
      const Complex V_2 = 8.*II*omega*rad + lambda;
      
      return(V_1 + V_2);
    };

  inline Complex Delta_times_V(const Real rad)
    {
      const Real Kay = K(rad);
      
      const Complex V_1 = -Kay*(Kay + 4.*II*(rad - 1.));
      const Complex V_2 = Kerr::Delta(rad, a)*(8.*II*omega*rad + lambda);
      
      return(V_1 + V_2);
    };
  
  inline Complex beta(Real const rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      
      return(2.*Delta*(-II*K(rad) + rad - 1. - 2.*Delta/rad));
    };

  inline Complex beta_over_Delta(Real const rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      
      return(2.*(-II*K(rad) + rad - 1. - 2.*Delta/rad));
    };

  inline Complex dr_beta(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      
      const Complex dr_beta_1 = 2.*(rad - 1. - 2.*Delta/rad -
				    II*Kay)*dr_Delta;
      const Complex dr_beta_2 = 2.*(1. + 2.*Delta/(rad*rad) -
				    2.*dr_Delta/rad - II*dr_Kay)*Delta;
      
      return(dr_beta_1 + dr_beta_2);
    };

  inline Complex ddr_beta(const Real rad)
    {
      const Real x = 1./rad;
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Real ddr_Delta = Kerr::ddr_Delta();
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real ddr_Kay = ddr_K();
      
      const Complex ddr_beta_1 = 4.*dr_Delta*(1. - 2.*x*(dr_Delta -
							 x*Delta) -
					      II*dr_Kay);
      const Complex ddr_beta_2 = 2.*(rad - 1. - 2.*x*Delta -
				     II*Kay)*ddr_Delta;
      const Complex ddr_beta_3 = 2.*Delta*(x*(-2*ddr_Delta +
					      x*(4.*dr_Delta +
						 x*(-4*Delta))) -
					   II*ddr_Kay);
      
      return(ddr_beta_1 + ddr_beta_2 + ddr_beta_3);
    };

  inline Complex alpha(const Real rad)
    {
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real Delta = Kerr::Delta(rad,a);
      const Complex bbeta = beta(rad);
      
      return(-II*Kay*bbeta/(Delta*Delta) + 3.*II*dr_Kay + lambda +
	     6.*Delta/(rad*rad));
    };

  Complex dr_alpha(const Real rad)
    {
      const Real x = 1./rad;
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real ddr_Kay = ddr_K();
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex bbeta = beta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      
      const Real dr_alpha_1 = 6.*x*x*(dr_Delta - 2.*x*Delta);
      const Complex dr_alpha_2 = -II*(Kay*dr_bbeta +
				      dr_Kay*bbeta)/(Delta*Delta);
      const Complex dr_alpha_3 = II*(2.*bbeta*Kay*dr_Delta/
				     (Delta*Delta*Delta) + 3.*ddr_Kay);
      
      return(dr_alpha_1 + dr_alpha_2 + dr_alpha_3);
    };
  
  Complex Deltasqr_times_dr_alpha(const Real rad)
    {
      const Real x = 1./rad;
      const Real Kay = K(rad);
      const Real dr_Kay = dr_K(rad);
      const Real ddr_Kay = ddr_K();
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex bbeta = beta(rad);
      const Complex bbeta_over_Delta = beta_over_Delta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      
      const Real dr_alpha_1 = 6.*x*x*(dr_Delta - 2.*x*Delta)*Delta*Delta;
      const Complex dr_alpha_2 = -II*(Kay*dr_bbeta + dr_Kay*bbeta);
      const Complex dr_alpha_3 = II*(2.*bbeta_over_Delta*Kay*dr_Delta +
				     3.*Delta*Delta*ddr_Kay);
      
      return(dr_alpha_1 + dr_alpha_2 + dr_alpha_3);
    };
  
  inline Complex U1(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex Vee = V(rad);
      const Complex Kay = K(rad);
      const Complex aalpha = alpha(rad);
      const Complex dr_aalpha = dr_alpha(rad);
      const Complex bbeta = beta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      const Complex ddr_bbeta = ddr_beta(rad);
      const Complex eeta = eta(rad);
      const Complex dr_eeta = dr_eta(rad);
      
      const Complex U1_1 = (2.*Delta*Delta*dr_aalpha -
			    dr_bbeta*dr_Delta)/bbeta;
      const Complex U1_2 = -Delta*dr_eeta*(aalpha*Delta +
					   dr_bbeta)/(bbeta*eeta);
      const Complex U1_3 = Delta*ddr_bbeta/bbeta;

      return(Vee + U1_1 + U1_2 + U1_3);
    };

  inline Complex Delta_times_U1(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real dr_Delta = Kerr::dr_Delta(rad);
      const Complex Vee = Delta_times_V(rad);
      const Complex Kay = K(rad);
      const Complex aalpha = alpha(rad);
      const Complex Deltasqr_times_dr_aalpha =
	Deltasqr_times_dr_alpha(rad);
      const Complex bbeta_over_Delta = beta_over_Delta(rad);
      const Complex dr_bbeta = dr_beta(rad);
      const Complex ddr_bbeta = ddr_beta(rad);
      const Complex eeta = eta(rad);
      const Complex dr_eeta = dr_eta(rad);
      
      const Complex U1_1 = (2.*Deltasqr_times_dr_aalpha -
 			    dr_bbeta*dr_Delta)/bbeta_over_Delta;
      const Complex U1_2 = -Delta*dr_eeta*(aalpha*Delta +
					   dr_bbeta)/(bbeta_over_Delta*eeta);
      const Complex U1_3 = Delta*ddr_bbeta/bbeta_over_Delta;

      return(Vee + U1_1 + U1_2 + U1_3);
    };

  inline Complex F(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real r2pa2 = rad*rad + a*a;
      const Complex eeta = eta(rad);
      const Complex dr_eeta = dr_eta(rad);
      
      const Complex tmp = (dr_eeta/eeta)*(Delta/r2pa2);
      
      return(tmp);
    };

  inline Complex U(const Real rad)
    {
      const Real Delta = Kerr::Delta(rad,a);
      const Real r2pa2 = rad*rad + a*a;
      const Real Gee = G(rad);
      const Real dr_Gee = dr_G(rad);
      const Complex DeltaYoo1 = Delta_times_U1(rad);
      const Complex Eff = F(rad);
      
      const Complex U_1 = DeltaYoo1/(r2pa2*r2pa2) + Gee*Gee;
      const Complex U_2 = Delta*dr_Gee/r2pa2 - Eff*Gee;
      
      return(U_1 + U_2);
    };

  Complex CalcR(const Complex X, const Complex dr_X, const Real rad);
  Complex Calcdr_R(const Complex X, const Complex dr_X,
		   const Complex ddr_X, const Real rad);
  Complex Calcddr_R(Complex R, Complex dr_R, const Real rad);

  Complex P(const Real rad);
  Complex dr_P(const Real rad);

  Real drsdr(const Real r);

  //
  // In what follows, Complex X[] has the following fields:
  // X[1].real() = r(star); X[1].imag() = 0.;
  // X[2] = Sasaki-Nakamura X;
  // X[3] = Sasaki-Nakamura dX/drstar;
  //
  // Real Y[] has the following fields:
  // Y[1] = r(rstar);
  // Y[2] = xi;
  // Y[3] = dxi/drstar;
  // Y[4] = psi;
  // Y[5] = dpsi/drstar;
  //
  // Sasaki-Nakamura X = exp(xi + II*psi);
  //
  void XInfInitCond(Complex XInf_Init[], const Real rout, const Real rsout);
  void XHInitCond(Real YH_Init[], const Real deltar, const Real rsin);
  void CalcXH(Real Y[], Real Yerr[], const Real rad);
  void CalcXInf(Complex XInf[], Complex XInferr[], const Real rad);
  Complex CalcAin(const Real rad);
  //
  // Overload the integration routines.
  //
  void IntegrateODE(Real y[], const Real rs1, const Real rs2,
		    const Real tol);
  void IntegrateODE(Complex y[], const Real rs1, const Real rs2,
		    const Real tol);
  void SNT_derivs(const Real rs, Real Y[], Real drs_Y[]);
  void SNT_derivs(const Real rs, Complex X[], Complex drs_X[]);
  //
  // Makes smooth phase functions from arctangents: ensures that
  // the answer doesn't deviate from a previous value by a factor of
  // 2 Pi.
  //
  Real my_atan2(const Real y, const Real x, const Real prev);
  //
  void odeint(Real ystart[], const int nvar, const Real x1,
	      const Real x2, const Real eps, const Real h1,
	      const Real hmin);
  void odeint(Complex ystart[], const int nvar, const Real x1,
	      const Real x2, const Real eps, const Real h1,
	      const Real hmin);
  void bsstep(Real y[], Real dydx[], const int nv,
	      Real *xx, Real htry, const Real eps,
	      Real yscal[], Real *hdid, Real *hnext);
  void bsstep(Complex y[], Complex dydx[], const int nv,
	      Real *xx, Real htry, const Real eps,
	      Real yscal[], Real *hdid, Real *hnext);
  void mmid(Real y[], Real dydx[], int nvar, Real xs,
	    Real htot, int nstep, Real yout[]);
  void mmid(Complex y[], Complex dydx[], int nvar, Real xs,
	    Real htot, int nstep, Complex yout[]);
};
#endif

