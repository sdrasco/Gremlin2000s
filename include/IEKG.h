//---------------------------------------------------------------------------
//
// $Id: IEKG.h,v 1.27 2007/07/18 22:27:54 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for inclined eccentric Kerr geodesics.
//
// [1] W. Schmidt Class. Quantum Grav. 19 2743 (2002) 
// [2] S. A. Hughes, Phys. Rev. D 61 084004 (2000)
//
// Steve Drasco
//
#ifndef _IEKG_H
#define _IEKG_H
#include <math.h>
#include "NRIEKG.h"

class IEKG {
public:
  // constructor
  IEKG(const int ProGrade, const Real spin, const Real eccentricity, 
       const Real SemiLatusRectum, const Real ThetaMinus, const Real tol);

  // overloaded constructor: this one lets you use E, Lz, and Q to specify the orbit
  IEKG(const Real spin, const Real E_in, const Real Lz_in, const Real Q_in, const Real tol);

  // destructor
  ~IEKG();

  // constants set by the constructor  
  int   prograde;         // 1 for prograde, 0 for retrograde
  Real	a;		  // spin parameter / M  ( 0 <=  a <= 1 )
  Real	e;		  // eccentricity
  Real	p;		  // semi latus rectum
  Real  theta_;           // theta turning point
  Real  eps;              // error tollerance
  Real	E;		  // Energy / mu
  Real	Lz;		  // z-component of angular momentum / (mu M)
  Real	Q;		  // Carter constant / (mu^2 M^2)
  Real	CosIota;	  // cos(iota) as defined in [2]
  Real	rmin;		  // inner radial turning point
  Real	rmax;		  // outter radial turning point
  Real	z_;		  // cos(theta_)
  Real  zPlusSquared;     // non-physical cos^2(theta) for theta_dot=0
  Real  zMinusSquared;    // cos^2(theta_)
  Real  ellK;             // common elliptic integral
  Real  beta;             // a*sqrt(1.0 - E*E)
  Real  betaRootzPlusSquared; // behaves badly in a=0 limit if not computed carefully (which it is)
  Real  aazPlusSquared;   // a^2 times zPlusSquared (same story as betaRootzPlusSquared)
  Real  zmzp;             // zMinusSquared / zPlusSquared (same story as betaRootzPlusSquared)
  Real  gamma;            // converts proper frequencies to BL frequencies
  Real	OmegaR;		  // M*(radial frequency in coordinate time)
  Real	OmegaTheta;	  // M*(theta frequency in coordniate time)
  Real	OmegaPhi; 	  // M*(phi frequency in coordinate time)
  Real  Gamma;            // Converts BL frequencies to Mino frequencies
  Real  UpsilonR;         // M*(radial frequency in Mino time)
  Real  UpsilonTheta;     // M*(theta frequency in Mino time)
  Real  UpsilonPhi;       // M*(Phi frequency in Mino time)
  int   theta_kmax;       // theta = theta_0 + 2 \sum{k=1}^theta_kmax ...
  int   r_nmax;           // theta = r_0 + 2 \sum{r=1}^r_nmax ...
  int   Delta_t_kmax;     // Delta_t = \sum_{k=1}^Delta_t_kmax ...
  int   Delta_t_nmax;     // Delta_t = ... \sum_{n=1}^Delta_t_nmax ...
  int   Delta_phi_kmax;   // same as two lines above with t -> phi
  int   Delta_phi_nmax;   // same as two lines above with t -> phi
  Real  *theta_k;         // array of 1 + theta_kmax coefficients
  Real  *r_n;             // array of 1 + r_nmax coefficients
  Real  *Delta_t_k;       // array of Delta_t_kmax coefficients
  Real  *Delta_t_n;       // array of Delta_t_nmax coefficients 
  Real  *Delta_phi_k;     // array of Delta_phi_kmax coefficients
  Real  *Delta_phi_n;     // array of Delta_phi_nmax coefficients
  int   k;                // index passed to integrands
  int   n;                // index passed to integrands
  Real  r_isco;           // radius of innermost stable circular orbit (ISCO)
  Real  Omega_isco;       // orbital frequency of ISCO
  Real  avgcos2;          // average of cos^2 theta
  Real  avgcot2;          // average of cot^2 theta
  
  // m*Omega_phi + k*Omega_theta + n*Omega_r
  Real  Omega(int mm, int kk, int nn); 

  // m*Upsilon_phi + k*Upsilon_theta + n*Upsilon_r
  Real  Upsilon(int mm, int kk, int nn); 

  // Delta_t(w_r,w_theta)
  Real  Delta_t(Real w_r, Real w_theta);
  Real  Delta_t(Real SgnRDot, Real SgnThetaDot, 
                Real w_r, Real w_theta);

  // Delta_phi(chi,psi)
  Real  Delta_phi(Real w_r, Real w_theta);
  Real  Delta_phi(Real SgnRDot, Real SgnThetaDot,
                  Real w_r, Real w_theta);

  // coordinates as a funciton of chi,psi
  Real Theta(Real chi);
  Real RofPsi(Real psi);

  // w_r, w_theta (and derivatives) as functions of chi and psi
  Real wTheta(Real chi);
  Real wR(Real psi);
  Real dwdchi(Real chi);
  Real dwdpsi(Real psi);

  // dt/dlambda as a (split) function of chi and psi
  Real TofChi(Real chi);
  Real TofPsi(Real psi);

private:

  // compute <cos^2 theta> and <cot^2 theta>
  void compute_avgcos2();
  Real ignd_avgcos2(Real chi);
  void compute_avgcot2();
  Real ignd_avgcot2(Real chi);

  // compute r_isco and Omega_isco 
  void compute_isco();

  // compute coordinate, proper, and Mino frequencies
  void  compute_freqs();
  
  // We will later evaluate many integrals of the form
  // int_0^x dx Xfunc(x).  It's easier to evaluate these
  // integrals using Chebyshev approximation of the integrand.
  void  prepare_Xfunc();

  // prepare coefficients for evaluating non-oscilating parts of t and phi (Deltat and Deltaphi)
  void  prepare_Deltas();

  // prepare coefficients for evaluating r(w_r) and theta(w_theta)
  void prepare_rtheta();

  // expansion coefficients
  Real  theta_kcoef(int kk);
  Real  r_ncoef(int nn);
  Real  Delta_t_kcoef(int kk);
  Real  Delta_t_ncoef(int nn);
  Real  Delta_phi_kcoef(int kk);
  Real  Delta_phi_ncoef(int nn);

  // integration routines
  Real qromb(Real (IEKG::*func)(Real),Real start, Real finish);
  void trpzd(Real (IEKG::*func)(Real),
              Real start,Real finish,int nn,Real last,Real *next);
  void polint(Real *xa,Real *ya,Real x,Real &y,Real &dy,int nn);
  Real qgaus(Real (IEKG::*func)(Real), Real start, Real finish);

  // dynamical step-size integration routines
  Real odeint(Real (IEKG::*derivs)(Real), Real x1, Real x2, Real h1);
  void rkqs(Real &y, Real dydx, Real &x, Real htry,
        Real eps, Real yscal, Real &hdid, Real &hnext,
        Real (IEKG::*derivs)(Real) );
  void rkck(Real &y, Real &dydx, Real x,
        Real h, Real &yout, Real &yerr,
        Real (IEKG::*derivs)(Real) );

  // Chebyshev Integration stuff
  Real *cint;
  int  CN;
  void chebft(Real aa, Real b, Real *c, int nn, Real (IEKG::*func)(Real));
  Real chebev(Real aa, Real b, Real *c, int m, Real x);
  void chint(Real aa, Real b, Real *c, Real *cint, int nn);

  Real F(Real x);
  Real G(Real x);
  Real H(Real x);
  Real J(Real x);
  Real Wfunc(Real x);
  Real Xfunc(Real x);
  Real Yfunc(Real x);
  Real Zfunc(Real x);
  Real PhiofChi(Real chi);
  Real PhiofPsi(Real psi);
  Real lambda0(Real chi);
  Real theta_k_integrand(Real chi);
  Real r_n_integrand(Real psi);
  Real Delta_t_k_integrand(Real chi);
  Real Delta_t_n_integrand(Real psi);
  Real Delta_phi_k_integrand(Real chi);
  Real Delta_phi_n_integrand(Real psi);

  // computes rmin, rmax, e, and p from a, E, Lz, Q.
  void e_and_p();
  void zroots(Complex a[], int m, Complex roots[], int polish, Real EPS);
  void laguer(Complex a[], int m, Complex *x, int *its);

};

#endif
