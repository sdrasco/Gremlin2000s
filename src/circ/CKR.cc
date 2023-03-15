//---------------------------------------------------------------------------
//
// $Id: CKR.cc,v 1.6 2004/09/16 21:51:35 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods of class CKR --- the class containing objects relevant
// to radiation for circular, non-equatorial orbits.
//
// Note that when Q==0 (equatorial orbits) the numerical behavior of
// the quantities that go into the definition of the frequencies
// \Omega_\theta etc is rather nasty: T_\theta diverges, \Omega_\phi
// is of the form infinity/infinity, no one is happy.  Q==0 should
// thus be treated as a special case.  This code does NOT do that;
// the CEKR class is designed for that case.
//
// NOT TRUE: those quantities are all well-haved at Q==0.  The
// problem is that the amplitude of motion goes to zero, and the
// parameters which are passed into the elliptic integrals were
// of the form 0/0.  By redefining these parameters, that problem
// no longer exists and this code can now duplicate the results
// of CEKR.cc.
//
// Inherits class Integral.
//
// Scott Hughes, culled from rad_circ.cc 9 January 1999.
//
/////////////////////////////////////////////////////////////////////
//
// MODIFIED: 5 April 1999
// I've changed the meaning of n, subdividing it into two different
// labels.  Motion for which thetadot is positive (ie from thetamin
// to thetamax) is seg == 0; thetadot negative is seg == 1.  In each
// segment, n == 0 denotes theta between thetamin and the equator,
// and n == 1 denotes theta between thetamax and the equator.
//
/////////////////////////////////////////////////////////////////////
//
// HEAVILY MODIFIED: 10 - 15 July 2004
// Integration methods have been largely reworked.  Much of the
// slowness in evaluating the integrals for ZH and ZInf comes from
// resolving wiggles in the integrand.  qromb() is crappy at this when
// there are many wiggles!
//
// New method: Expand the integrand in Chebyshev polynomials, and then
// calculate the integral by Clenshaw-Curtis quadrature.  This amounts
// to a simple sum over the coefficients of the expansion --- see
// Numerical Recipes (2nd edition), end of Section 5.9, esp.  Equation
// (5.9.3).  Further, take advantage of the fact that the integral we
// need to evaluate the expansion coefficients is just a cosine
// Fourier transform --- see Numerical Recipes (2nd edition), Section
// 5.8, esp. text following the "chebft()" code.  Use cosft2() for
// this transform.
//
// Numerical experiments show that this speeds up the code by a large
// factor --- we don't hang up on evaluation of Zeds with wiggly
// integrands anymore.
//
// Code no longer inherits Integral.
//
// MINOR MODIFICATIONS: 16 - 17 July 2004
// Made the integrator adaptive.  This allows us to estimate the
// accuracy of integration (not possible before); and, it improves the
// behavior when (l,m,k) are all largish.  The "vanilla" C-C
// quadrature was exhibiting crappy behavior for some orbits with, for
// example, l = 9, m = 9, k > 15.  Improving this was simple --- just
// include more coefficients in the expansion --- but very ineffecient
// for the cases that were already working well.  The new algorithm,
// which adapative increases ncof and stops when the results are
// "accurate enough" (or, until we've iterated too many times; in
// those cases, we typically have very small values of ZedH, ZedInf).
//
// In order to maximally re-use previous calculations, we now use the
// "trapezoidal" Gauss-Lobatto rule for the expansion coefficients,
// cf. NumRec ed 2 Eq (5.9.4).  (Be careful with indexing!)  This
// means that we use cosft1() for our cosine transform, which in turn
// means that our initial storage is for ncofs + 1 points rather than
// ncofs.  The C-C quadrature is otherwise identical.  The main new
// element is that we perform the quadrature at least twice, for ncofs
// = 64 and ncofs = 128.  We then use ratint() to extrapolate to ncofs
// = infinity (1/ncofs = 0).
//
#include <math.h>
#include "Globals.h"
#include "CKR.h"
#include "SWSH.h"
#include "CKG.h"
#include "SNT.h"

// Constructor.
CKR::CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon) :
  swsh(swsh_in), snt(snt_in), ckg(ckg_in)
{
  //
  // Local storage of orbital characteristics.
  //
  r = ckg->r;
  a = ckg->a;
  E = ckg->E;
  Lz = ckg->Lz;
  Q = ckg->Q;
  //
  // Local storage of mode characteristics.
  //
  l = swsh->Gl;
  m = ckg->m;
  k = ckg->k;
  wmk = ckg->wmk;
  pmk = ckg->pmk;
  //
  // Local storage of useful quantities that appear all over the
  // place.
  //
  Delta = Kerr::Delta(r, a);
  dr_Delta = Kerr::dr_Delta(r);
  SN_Kay = snt->K(r);
  dr_SN_Kay = snt->dr_K(r);
  //
  ComputeZeds(ZedH, ZedInf, epsilon);
}

// The phase.
Real CKR::Psi_mk(const Real chi)
{
  return(wmk*ckg->t(chi) - m*ckg->phi(chi));
}

// The Newman-Penrose quantity rho (with sign as in
// Teukolsky's papers).
Complex CKR::rho_func(const int n, const Real z)
{
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Complex tmp = Complex(r, -a*ct);

  return(-1./tmp);
}

// The C_{ab} functions are projections of the orbiting
// particle's stress-energy tensor onto the Newman-Penrose
// tetrad (modulo some factors).
Complex CKR::Cnn(const int seg, const int n, const Real z)
{
  const Real tmp = E*(r*r + a*a) - a*Lz;
  const Real sig = Kerr::Sigma(r, a, z);

  return(0.25*tmp*tmp/(sig*sig*gkg.TFunc(r, a, z, E, Lz, Q)));
}

Complex CKR::Cnmbar(const int seg, const int n, const Real z)
{
  const Real st = sqrt(1.-z);
  const Complex rho = rho_func(n, z);
  const Real sig = Kerr::Sigma(r, a, z);

  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  Complex tmp2;

  if(seg == 0)
    tmp2 = Complex(sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), st*a*E - Lz/st);
  else
    tmp2 = Complex(-sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), st*a*E - Lz/st);

  return(rho*tmp1*tmp2/(2.*sqrt(2.)*sig*gkg.TFunc(r, a, z, E, Lz, Q)));
}

Complex CKR::Cmbarmbar(const int seg, const int n, const Real z)
{
  const Real st = sqrt(1. - z);
  const Complex rho = rho_func(n, z);

  Complex tmp;

  if (seg == 0)
    tmp = Complex(sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), a*E*st - Lz/st);
  else
    tmp = Complex(-sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), a*E*st - Lz/st);

  return(rho*rho*tmp*tmp/(2.*gkg.TFunc(r, a, z, E, Lz, Q)));
}

// The A_{abj}'s are functions that appear when the source is
// integrated over the Green's function.  See notes for details
// --- it's a touch involved.
Complex CKR::Ann0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);

  const Complex tmp1 = 2.*II*a*rho*st*swsh->l2dagspheroid(ct);
  const Real tmp2 = swsh->l1dagl2dagspheroid(ct);

  const Complex pref = -2.*Cnn(seg, n, z)/(Delta*Delta*rhob*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CKR::Anmbar0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);
  const Real S = swsh->spheroid(ct);
  const Real L2dagS = swsh->l2dagspheroid(ct);

  const Complex tmp1 = (II*SN_Kay/Delta - rho - rhob)*L2dagS;
  const Complex tmp2 = (II*SN_Kay/Delta + rho +
			rhob)*II*a*st*S*(rho - rhob);

  const Complex pref = -2.*sqrt(2.)*Cnmbar(seg, n, z)/(Delta*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CKR::Anmbar1(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);
  const Real S = swsh->spheroid(ct);
  const Real L2dagS = swsh->l2dagspheroid(ct);

  const Complex tmp = L2dagS + II*a*st*(rho - rhob)*S;

  const Complex pref = -2.*sqrt(2.)*Cnmbar(seg, n, z)/(Delta*rho3);
  return(pref*tmp);
}

Complex CKR::Ambarmbar0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  const Complex tmp1 = (SN_Kay/Delta)*(SN_Kay/Delta);
  const Complex tmp2 = 2.*rho*II*SN_Kay/Delta;
  const Complex tmp3 = II*(dr_SN_Kay - SN_Kay*dr_Delta/Delta)/Delta;

  const Complex pref = S*Cmbarmbar(seg, n, z)*rhob/rho3;

  return(pref*(tmp1 + tmp2 + tmp3));
}

Complex CKR::Ambarmbar1(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  const Complex tmp = -II*SN_Kay/Delta;

  const Complex pref = 2.*S*rhob*Cmbarmbar(seg, n, z)/rho3;

  return(pref*(rho + tmp));
}

Complex CKR::Ambarmbar2(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  return(-Cmbarmbar(seg, n, z)*rhob*S/rho3);
}

// "dZed" is what you get when you do the radial integral of
// source and Green's function.  Integrate it and you get
// Zed [modulo factors which are handled in the function
// Zed().]
Complex CKR::dZedH(const int seg, const int n, const Real chi)
{
  const Real cchi = cos(chi);
  const Real z = ckg->zedminus*cchi*cchi;
  const Real sign = (seg == 0) ? 1. : -1.;

  const Complex term0 = snt->TeukRH*(Ann0(seg, n, z) +
				     Anmbar0(seg, n, z) +
				     Ambarmbar0(seg, n, z));
  const Complex term1 = -snt->dr_TeukRH*(Anmbar1(seg, n, z) +
					 Ambarmbar1(seg, n, z));
  const Complex term2 = snt->ddr_TeukRH*Ambarmbar2(seg, n, z);

  const Complex pref = ckg->dtdchi(chi)*exp(sign*II*Psi_mk(chi));

  return(pref*(term0 + term1 + term2));
}

// A wrapper for dZedH.
Complex CKR::IntegrandH(const Real chi)
{
  const int n = (chi < 0.5*M_PI ? 0 : 1);

  return(dZedH(0, n, chi) + dZedH(1, n, chi));
}

Complex CKR::dZedInf(const int seg, const int n, const Real chi)
{
  const Real cchi = cos(chi);
  const Real z = ckg->zedminus*cchi*cchi;
  const Real sign = (seg == 0) ? 1. : -1.;

  const Complex term0 = snt->TeukRInf*(Ann0(seg, n, z) +
				       Anmbar0(seg, n, z) +
				       Ambarmbar0(seg, n, z));
  const Complex term1 = -snt->dr_TeukRInf*(Anmbar1(seg, n, z) +
					   Ambarmbar1(seg, n, z));
  const Complex term2 = snt->ddr_TeukRInf*Ambarmbar2(seg, n, z);

  const Complex pref = ckg->dtdchi(chi)*exp(sign*II*Psi_mk(chi));

  return(pref*(term0 + term1 + term2));
}

// A wrapper for dZedInf.
Complex CKR::IntegrandInf(const Real chi)
{
  const int n = (chi < 0.5*M_PI ? 0 : 1);

  return(dZedInf(0, n, chi) + dZedInf(1, n, chi));
}

void CKR::ComputeZeds(Complex & ZedH, Complex & ZedInf, const Real epsilon)
{
  //
  // Compute amplitudes ZedH, ZedInf ... Uses new algorithm based on
  // Chebyshev expansion followed by Clenshaw-Curtis quadrature.
  // Further modified to be adaptive in hopes that this will fix up
  // some of the convergence problems I've been seeing with l, m, and
  // k are all large.
  //
  int iter, itmax = 15;
  int i, ncofs = 64;
  Complex Hsum, Infsum;
  Real *oldy, *newy;
  Real pio2 = 0.5*M_PI;
  Complex *oldHcofs, *newHcofs;
  Complex *oldInfcofs, *newInfcofs;
  Complex maxHcof, maxInfcof;
  //
  // Storage and utility variables for the rational extrapolation.
  //
  int K = 2;
  Complex HSum, InfSum, HSumerr, InfSumerr;
  Complex *Hsumest;    Hsumest = Complexvector(1, itmax);
  Complex *Infsumest;  Infsumest = Complexvector(1, itmax);
  Real *oneovern;      oneovern = Realvector(1, itmax);
  for (iter = 1; iter <= itmax; iter++) {
    if (iter == 1) {
      newy = Realvector(1, ncofs + 1);
      oldy = Realvector(1, ncofs + 1);
      newHcofs = Complexvector(1, ncofs + 1);
      oldHcofs = Complexvector(1, ncofs + 1);
      newInfcofs = Complexvector(1, ncofs + 1);
      oldInfcofs = Complexvector(1, ncofs + 1);
      for (i = 1; i <= ncofs + 1; i++) {
	//
	// Evaluating the integrand at the proper points in prep for
	// the Chebyshev expansion.  See NR (2nd ed), Eqs (5.9.4) and
	// (5.8.10) (with a = 0, b = pi/2).
	//
	newy[i] = cos(M_PI*(i - 1)/((Real)ncofs));
	newHcofs[i] = IntegrandH(newy[i]*pio2 + pio2);
	newInfcofs[i] = IntegrandInf(newy[i]*pio2 + pio2);
	//
	// Copy data into "old" arrays to be used in the next
	// iteration.
	//
	oldy[i] = newy[i];
	oldHcofs[i] = newHcofs[i];
	oldInfcofs[i] = newInfcofs[i];
      }
      //
      // The cosine Fourier transform.  Note, I have modified the NR
      // FFT routines very slightly to handle complex functions.
      // Since there is no crosstalk between the real and the
      // imaginary parts, this is equivalent to FFT'ing a pair of real
      // functions.  I found experimentally that doing it this way is
      // quite a bit faster than storing the real and imaginary parts
      // separately and performing two FFTs.
      //
      cosft1(newHcofs, ncofs);
      cosft1(newInfcofs, ncofs);
      //
      // Clenshaw-Curtis quadrature --- NR ed 2, Eq (5.9.3).
      maxHcof = newHcofs[1];
      maxInfcof = newInfcofs[1];
      Hsum = 0.5*newHcofs[1];
      Infsum = 0.5*newInfcofs[1];
      for (i = 1; i <= ncofs/2; i++) {
	Hsum -= newHcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
	Infsum -= newInfcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
	//
	if (abs(newHcofs[2*i + 1]) > abs(maxHcof))
	  maxHcof = newHcofs[2*i + 1];
	if (abs(newInfcofs[2*i + 1]) > abs(maxInfcof))
	  maxInfcof = newInfcofs[2*i + 1];
      }
    } else {
      //
      // This bracketed section is identical to the above section,
      // except that we carefully reuse data that was calculated on
      // the previous iteration whenever possible.
      //
      ncofs *= 2;
      newy = Realvector(1, ncofs + 1);
      newHcofs = Complexvector(1, ncofs + 1);
      newInfcofs = Complexvector(1, ncofs + 1);
      //
      // Points we can recycle from the previous iteration
      for (i = 1; i <= ncofs + 1; i += 2) {
	newy[i] = oldy[(i - 1)/2 + 1];
	newHcofs[i] = oldHcofs[(i - 1)/2 + 1];
	newInfcofs[i] = oldInfcofs[(i - 1)/2 + 1];
      }
      //
      // Points we cannot recycle
      for (i = 2; i <= ncofs; i += 2) {
	newy[i] = cos(M_PI*(i - 1)/((Real)ncofs));
	newHcofs[i] = IntegrandH(newy[i]*pio2 + pio2);
	newInfcofs[i] = IntegrandInf(newy[i]*pio2 + pio2);
      }
      //
      // Remake the "old" arrays
      free_Realvector(oldy, 1, ncofs/2 + 1);
      oldy = Realvector(1, ncofs + 1);
      free_Complexvector(oldHcofs, 1, ncofs/2 + 1);
      oldHcofs = Complexvector(1, ncofs + 1);
      free_Complexvector(oldInfcofs, 1, ncofs/2 + 1);
      oldInfcofs = Complexvector(1, ncofs + 1);
      for (i = 1; i <= ncofs + 1; i++) {
	oldy[i] = newy[i];
	oldHcofs[i] = newHcofs[i];
	oldInfcofs[i] = newInfcofs[i];
      }
      cosft1(newHcofs, ncofs);
      cosft1(newInfcofs, ncofs);
      //
      maxHcof = newHcofs[1];
      maxInfcof = newInfcofs[1];
      Hsum = 0.5*newHcofs[1];
      Infsum = 0.5*newInfcofs[1];
      for (i = 1; i <= ncofs/2; i++) {
	Hsum -= newHcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
	Infsum -= newInfcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
	if (abs(newHcofs[2*i + 1]) > abs(maxHcof))
	  maxHcof = newHcofs[2*i + 1];
	if (abs(newInfcofs[2*i + 1]) > abs(maxInfcof))
	  maxInfcof = newInfcofs[2*i + 1];
      }
    }
    Hsum *= 2.*M_PI/((Real)ncofs);
    Infsum *= 2.*M_PI/((Real)ncofs);
    maxHcof *= 2.*M_PI/((Real)ncofs);
    maxInfcof *= 2.*M_PI/((Real)ncofs);
    //
    free_Realvector(newy, 1, ncofs + 1);
    free_Complexvector(newHcofs, 1, ncofs + 1);
    free_Complexvector(newInfcofs, 1, ncofs + 1);
    //
    if (abs(Hsum)/abs(maxHcof) < 1.e-11 &&
	abs(Infsum)/abs(maxInfcof) < 1.e-11) {
      //
      // The sum is so close to zero we're just getting roundoff crap.
      // Set it to zero and break.
      Hsum = 0.; Infsum = 0.; break;
    }
    Hsumest[iter] = Hsum;
    Infsumest[iter] = Infsum;
    oneovern[iter] = 1./((Real)ncofs);
    if (iter >= K) {
      ratint(&oneovern[iter - K], &Hsumest[iter - K], K, 0., &HSum,
 	     &HSumerr);
      ratint(&oneovern[iter - K], &Infsumest[iter - K], K, 0., &InfSum,
 	     &InfSumerr);
      if (fabs(HSumerr.real()) <= epsilon*fabs(HSum.real()) &&
	  fabs(HSumerr.imag()) <= epsilon*fabs(HSum.imag()) &&
	  fabs(HSumerr.real()) <= epsilon*fabs(HSum.real()) &&
	  fabs(HSumerr.imag()) <= epsilon*fabs(HSum.imag())) break;
    }
  }
  // start test
  //
  // write a warning if we needed more terms.  (output as Matlab comment)
  /*
  if(iter >= itmax){
    cerr << "% WARNING CKR::ClenshawCurtisInt:  iter >= itmax = " << itmax
         << ".  We need more than " << ncofs << " terms.\n"
         << "% abs(sum)/abs(maxcof) = " << abs(Infsum)/abs(maxInfcof)
         << ";  sum = " << Infsum << endl;
    //for(int jj=1; jj<ncofs; jj++) cout << jj << "\t" << Infsumest[jj].real() << "\t" << sumest[jj].imag() << endl;
    //for(int jj=1; jj<ncofs; jj++) cout << jj << "\t" << oldInfcofs[jj].real() << "\t" << oldInfcofs[jj].imag() << endl;
    //Complex zz;
    //for(Real xx = 0.01; xx < M_PI; xx+= M_PI/100.0){
    //  zz = IntegrandInf(xx);
    //  cout << xx << "\t" << zz.real() << "\t" << zz.imag() << endl;
    //}
    exit(0);
  }
  */
  // end test
  const Complex denom = II*wmk*snt->Bin*ckg->T_th;
  ZedH = M_PI*Hsum/denom;
  ZedInf = M_PI*Infsum/denom;
  //
  free_Realvector(oldy, 1, ncofs + 1);
  free_Complexvector(oldHcofs, 1, ncofs + 1);
  free_Complexvector(oldInfcofs, 1, ncofs + 1);
  //
  free_Complexvector(Hsumest, 1, itmax);
  free_Complexvector(Infsumest, 1, itmax);
  free_Realvector(oneovern, 1, itmax);
}

// Cosine FFT, modified slightly for complex arguments.
void CKR::cosft1(Complex y[], int n)
{
  int j, n2;
  Complex sum, y1, y2;
  Real theta, wi=0.0, wpi, wpr, wr=1.0, wtemp;
  
  theta = M_PI/n;
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  sum = 0.5*(y[1] - y[n + 1]);
  y[1] = 0.5*(y[1] + y[n + 1]);
  n2 = n + 2;
  for (j = 2; j <= (n>>1); j++) {
    wr = (wtemp = wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
    y1 = 0.5*(y[j] + y[n2 - j]);
    y2 = (y[j] - y[n2 - j]);
    y[j] = y1 - wi*y2;
    y[n2 - j] = y1 + wi*y2;
    sum += wr*y2;
  }
  realft(y, n, 1);
  y[n + 1] = y[2];
  y[2] = sum;
  for (j = 4; j <= n; j += 2) {
    sum += y[j];
    y[j] = sum;
  }
}

// Real FFT, modified slightly for complex arguments.  (I guess it's
// not Real FFT anymore, eh?)
void CKR::realft(Complex data[], unsigned long n, int isign)
{
  unsigned long i, i1, i2, i3, i4, np3;
  Real c1 = 0.5, c2;
  Complex h1r, h1i, h2r, h2i;
  Real wr, wi, wpr, wpi, wtemp, theta;
  
  theta = 3.141592653589793/(Real) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, n>>1, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n+3;
  for (i = 2; i <= (n>>2); i++) {
    i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
    h1r = c1*(data[i1] + data[i3]);
    h1i = c1*(data[i2] - data[i4]);
    h2r = -c2*(data[i2] + data[i4]);
    h2i = c2*(data[i1] - data[i3]);
    data[i1] = h1r + wr*h2r - wi*h2i;
    data[i2] = h1i + wr*h2i + wi*h2r;
    data[i3] = h1r - wr*h2r + wi*h2i;
    data[i4] = -h1i + wr*h2i + wi*h2r;
    wr = (wtemp = wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
  }
  if (isign == 1) {
    data[1] = (h1r = data[1]) + data[2];
    data[2] = h1r - data[2];
  } else {
    data[1] = c1*((h1r = data[1]) + data[2]);
    data[2] = c1*(h1r - data[2]);
    four1(data, n>>1, -1);
  }
}

// Generic FFT, modified slightly for complex arguments.
void CKR::four1(Complex data[], unsigned long nn, int isign)
{
  unsigned long n, mmax, m, j, istep, i;
  Real wtemp, wr, wpr, wpi, wi, theta;
  Complex tempr, tempi;
  
  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i+= 2) {
    if (j > i) {
      tempr = data[j];
      data[j] = data[i];
      data[i] = tempr;
      //
      tempr = data[j + 1];
      data[j + 1] = data[i + 1];
      data[i + 1] = tempr;
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(2.*M_PI/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m+= 2) {
      for (i = m; i <= n; i+= istep) {
	j = i + mmax;
	tempr = wr*data[j] - wi*data[j + 1];
	tempi = wr*data[j + 1] + wi*data[j];
	data[j] = data[i] - tempr;
	data[j+1] = data[i + 1] - tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}

void CKR::ratint(const Real xa[], const Complex ya[], const int n,
		 const Real x, Complex *y, Complex *dy)
{
  int m, i, ns = 1;
  const Real TINY = 1.0e-25;
  Real hh, h;
  Complex w, t, dd, *c, *d;

  c = Complexvector(1, n);
  d = Complexvector(1, n);
  hh = fabs(x - xa[1]);
  for (i = 1; i <= n; i++) {
    h = fabs(x - xa[i]);
    if (h == 0.0) {
      *y = ya[i];
      *dy = 0.0;
      free_Complexvector(d, 1, n);
      free_Complexvector(c, 1, n);
      return;
    } else if (h < hh) {
      ns = i;
      hh = h;
    }
    c[i] = ya[i];
    d[i] = ya[i] + TINY;
  }
  *y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      w = c[i+1] - d[i];
      h = xa[i+m] - x;
      t = (xa[i] - x)*d[i]/h;
      dd = t - c[i+1] + TINY;
      if (dd.real() == 0.0 && dd.imag() == 0.0)
	Die("Error in routine ratint");
      dd = w/dd;
      d[i] = c[i+1]*dd + TINY;
      c[i] = t*dd;
    }
    *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_Complexvector(d, 1, n);
  free_Complexvector(c, 1, n);
  return;
}
