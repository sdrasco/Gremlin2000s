//
// Class IEKR (Inclined Eccentric Kerr Radiation).
//
// Steve Drasco
//
#ifndef _IEKR_H
#define _IEKR_H

#include "Globals.h"
#include "SWSH.h"
#include "GKG.h"
#include "IEKG.h"
#include "SNT.h"

class IEKR{
public:
  IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index, 
      int n_index, const Real tol);
  IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index,
      int n_index, const Real tol, Real MaxEIin, Real MaxEHin);
  IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index,
      int n_index, const Real tol, Real MaxEIin, Real MaxEHin, int DoZH, int DoZI);

  IEKG  *iekg; // Inclined Eccentric Kerr Geodesic object
  int       l; // Spheroidal harmonic l-index
  int       m; // Spheroidal harmonic m-index (orbital harmonic azimuthal index)
  int       k; // orbital harmonic polar index
  int       n; // orbital harmonic radial index
  Complex   ZH;// Z^H_{lmkn}  at r_f = \infty  (needed for fluxes at infinity)
  Complex   ZInf; // Z^{\infty}_{lmkn} at r_f = [horizon] needed for horizon fluxes
  Complex   Ain; // snt->Ain
  Real  omega; // Omega_{mkn}
  Real  alpha; // factor that appears in the fluxes on the horizon

private:

  SNT *snt;           // SNT object 
  SWSH *swsh;         // spin weighted spheroidal harmonic object
  Real a;             // spin parameter / M
  Real E;             // Energy / mu
  Real Lz;            // z-component of angular momentum / (mu M)
  Real Q;             // Carter constant / (mu^2 M^2)  
  Real Gamma;         // Converts BL frequencies to Mino frequencies
  Complex Almkn;      // A_{lmkn} in definition of Z^H_{lmkn} [see note in code]
  Complex AlmknInf;   // A_{lmkn} in definition of Z^\infty_{lmkn} [""]
  Real eps;           // error tollerance
  Real MaxEI;         // threshold value.  for modes with EI < (1e-9)*maxEI, ZH is good only to one digit
  Real MaxEH;         // threshold value.  for modes with EH < (1e-9)*maxEH, ZInf is good only to one digit
  GKG gkg;            // general kerr geodesic functions [R(r) and Theta(theta) etc]
  int IntToPi;        // 1 if integrals go from 0 to pi.  Otherwise integrals go from 0 to 2pi.
  int DoH;            // 1 if we are in the process of computing ZH, otherwise we do ZInf
  int psi_MinIterate; // minimum number of iterations for Romberg psi-integral 
  int chi_MinIterate; // minimum number of iterations for Romberg chi-integral
  Real psi_h_guess;   // initial stepsize for adaptive stepsize psi-integral
  Real chi_h_guess;   // initial stepsize for adaptive stepsize chi-integral
  int  ScaleIntegrand;// a switch to turn on scaling of the radial (psi) integrand
  Real ScaleLambda;   // the scaling is exp(-psi/lambda)

  // There are many parameters, such as say r(psi) which are used
  // during the integration by several different routines.  We could
  // write many short and sweet routines to compute these simple
  // quantitites in an effor to make the code for something like 
  // the A_{abc}'s very transperent.  Doing this however has the 
  // down side that the simple functions are called many times, and 
  // this increases the time it takes to do the full calculation.
  // So instead of simple routines, we will compute these quantities
  // inside the integrands themselves, oncer per "time step".  
  // The results  will be stored so that they can be used by, say, 
  // the A_{abc} routines.  Here we list those quantites:
  Real psi;         // integration varriable
  Real r;           // r(psi) 
  Real SgnRDot;     // 1 for psi < pi, -1 for psi > pi
  Real RDot;        // dr/dlambda, with sign above
  Real Delta;       // r^2 - 2r +a^2
  Real Sigma;       // r^2 + a^2 cos^2 theta
  Real KoverDelta;  // Delta / [(r^2 + a^2)omega - ma]
  Real dKoverDelta; // d/dr  (K/Delta)
  Real chi;         // integration varriable
  Real theta;       // theta(chi) 
  Real SgnThetaDot; // 1 for chi < pi, -1 for chi > pi
  Real ThetaDot;    // dtheta/dlambda, with sign above
  Real TDot;        // dt/dlambda
  Real ct;          // cos(theta) 
  Real st;          // sin(theta)
  Real cott;        // cot(theta) = cos(theta)/sin(theta)
  Real csct;        // csc(theta) = 1/sin(theta)
  Complex rho;      // -[r - ia cos(theta)]^{-1}
  Complex rhobar;   // -[r + ia cos(theta)]^{-1}
  Real wR;          // w_r(psi)
  Real wTheta;      // w_theta(chi)

  // items from the spin weighted spheroidal harmonic object 
  Complex S;     // S_{lmkn}(\theta)
  Complex L2S;   // L_2^\dagger S_{lmkn}(\theta)
  Complex L1L2S; // L_1^\dagger L_2^\dagger S_{lmkn}(\theta)

  // items from the Sasaki-Nakamura-Teukolsky object 
  // (for the homogeneous solution)
  Complex R;   // R^H_{lmkn}(r)
  Complex dR;  // d/dr R^H_{lmkn}(r)
  Complex ddR; // d/dr d/dr R^H_{lmkn}(r)
  Complex Bin;  // B_{in} 

  // Z_{lmnk} = (2\pi\Gamma)^{-1} \int_0^{2\pi} d\psi OuterPsiIntegrand
  Complex OuterPsiIntegrand(Real psi_tmp);

  // InnerPsiIntegrand = dw_r/d\psi(\psi) \exp^{i n w_r(\psi)} 
  //                     \times \int_0^{2\pi} d\chi ChiIntegrand
  Complex InnerChiIntegrand(Real chi_tmp);

  // Z_{lmnk} = (2\pi\Gamma)^{-1} \int_0^{2\pi} d\psi  \int_0^{2\pi} d\chi RawIntegrand(chi,psi)
  Complex RawIntegrand(Real chi_tmp, Real psi_tmp);

  // determines rescaling parameters: ScaleLambda_here = Psi_paper
  void EstimateScaleLambda();

  // subroutines contaning the innermost parts of the double integral
  Complex Jlm();
  Complex Ilm();

  // the A_{abc}'s
  Complex A_n_n_0();
  Complex A_n_mbar_0();
  Complex A_mbar_mbar_0();
  Complex A_n_mbar_1();
  Complex A_mbar_mbar_1();
  Complex A_mbar_mbar_2();

  // the C_{ab}'s
  Real C_n_n();
  Complex C_mbar_mbar();
  Complex C_n_mbar();

  // integration wrapper function
  Complex Integrate(Complex (IEKR::*igrnd)(Real), Real start, Real finish, int PsiIntegral, Real eps_int);

  // Romberg integration routines
  Complex RombergInt(Complex (IEKR::*func)(Real),Real start, Real finish, int MinIterate, Real eps_int);
  void Ctrpzd(Complex (IEKR::*func)(Real),
              Real start,Real finish,int nn,Complex last,Complex *next);
  void Cpolint(Real *xa,Complex *ya,Real x,Complex &y,Complex &dy,int nn);

  // Adaptive stepsize integration routines
  Complex CashCarpInt(Complex (IEKR::*derivs)(Real), Real x1, Real x2, Real h_guess, Real eps_int);
  void    Crkqs(Complex &y, Complex dydx, Real &x, Real htry, Real eps_int, 
            Real yscal, Real &hdid, Real &hnext, Complex (IEKR::*derivs)(Real) );
  void    Crkck(Complex &y, Complex &dydx, Real x, Real h, Complex &yout, 
            Complex &yerr, Complex (IEKR::*derivs)(Real) );

  // Clenshaw-Curtis integration routines
  Complex ClenshawCurtisInt(Complex (IEKR::*)(Real), Real x1, Real x2, int ncofs, int itmax, Real eps_int);
  void cosft1(Complex y[], int nn);
  void realft(Complex data[], unsigned long nn, int isign);
  void four1(Complex data[], unsigned long nn, int isign);
  void ratint(const Real xa[], const Complex ya[], const int n, const Real x, Complex *y, Complex *dy);

  // error function: erf(x) = (2/pi) int_0^x dy exp(-y^2)
  Real erff(Real x);
  Real gammp(Real aa, Real x);
  void gser(Real &gamser, Real aa, Real x, Real &gln);
  void gcf(Real &gammcf, Real aa, Real x, Real &gln);
  Real gammln(Real xx);

};

#endif
