//---------------------------------------------------------------------------
//
// $Id: CKG.cc,v 1.3 2004/09/17 19:33:05 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular Kerr geodesics; culled from geodesic.cc
// by Scott Hughes, 10 January 1999.
//
// MODIFIED: 5 April 1999
// I've changed the meaning of n, subdividing it into two different
// labels.  Motion for which thetadot is positive (ie from thetamin
// to thetamax) is seg == 0; thetadot negative is seg == 1.  In each
// segment, n == 0 denotes theta between thetamin and the equator,
// and n == 1 denotes theta between the thetamax and the equator.
//
#include <math.h>
#include "Globals.h"
#include "CKG.h"
#include "NRCKG.h"
#include "NRUtil.h"

#define EPSILON 1.e-9

CKG::CKG(const int emm, const int kay, const Real Lz_or_E_or_cosiota,
	 int Lz_or_E_or_cosiota_flag, const Real rad, const Real spin) :
  m(emm), k(kay), r(rad), a(spin)
{
  // This is needed right away ...
  Delta = Kerr::Delta(r, a);

  if (fabs(a) > 1e-7) {
    if (Lz_or_E_or_cosiota_flag == USING_Lz) {
      Lz = Lz_or_E_or_cosiota;
      E = Efunc();
      Q = Qfunc();
      cosiota = Lz/sqrt(Lz*Lz + Q);
    } else if (Lz_or_E_or_cosiota_flag == USING_E) {
      E = Lz_or_E_or_cosiota;
      Lz = Lzfunc();
      Q = Qfunc();
      cosiota = Lz/sqrt(Lz*Lz + Q);
    } else {
      cosiota = Lz_or_E_or_cosiota;
      //
      // Now, need to bisect this puppy until I get Lz and Q correct.
      //
      Real Ehi, Elo, cosiota_guess;
      Elo = Kerr::Eeqpro(r, a);
      if (r >= Kerr::isco_ret(a))
	Ehi = Kerr::Eeqret(r, a);
      else
	Ehi = 1.;
      int BRACKETED = 0, count = 0;
      while (!BRACKETED) {
	E = 0.5*(Ehi + Elo);
	Lz = Lzfunc();
	Q = Qfunc();
	cosiota_guess = Lz/sqrt(Lz*Lz + Q);
	if (cosiota_guess > cosiota) // Not energetic enough
	  Elo = E;
	else // too energetic
	  Ehi = E;
	count++;
	if (fabs(cosiota_guess - cosiota) < EPSILON) BRACKETED = 1;
 	else if (count > 200)
 	  Die("Unable to converge on stable orbit parameters!");
      }
      Lz = Lzfunc();
      Q = Qfunc();
      if (fabs(cosiota) < EPSILON) Lz = 0.;
    }
  } else {
    //
    // The above code breaks for Schwarzschild ... fortunately, that
    // case is quite simple!
    //
    const Real Ltot = Kerr::Lzeqpro(r, a);
    if (Lz_or_E_or_cosiota_flag == USING_Lz) {
      Lz = Lz_or_E_or_cosiota;
      E = Efunc();
      Q = Ltot*Ltot - Lz*Lz;
      cosiota = Lz/sqrt(Lz*Lz + Q);
    } else if (Lz_or_E_or_cosiota_flag == USING_E) {
      Die("Using E to specify orbital inclination for a = 0 doesn't work!");
    } else {
      cosiota = Lz_or_E_or_cosiota;
      Lz = Ltot*cosiota;
      Q = Ltot*Ltot - Lz*Lz;
      E = Efunc();
    }
  }
  alpha_wilkins = Q + Lz*Lz;
  beta_wilkins = a*a*(1. - E)*(1. + E);
  const Real tmp = alpha_wilkins + beta_wilkins;

  if (beta_wilkins > EPSILON) {
    zedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins));
    zedplus /= 2.*beta_wilkins;
    betazedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins))/2.;
    zedminus = (tmp - sqrt(tmp*tmp - 4.*Q*beta_wilkins));
    zedminus /= 2.*beta_wilkins;
    if (fabs(zedminus) < 1.e-14) zedminus = 0.0;
    k_zmzp = sqrt(zedminus/zedplus);
  } else {
    zedplus = tmp/beta_wilkins - (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    betazedplus = tmp - (Q*beta_wilkins/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    zedminus = (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    k_zmzp = 0.;
  }

  thetamin = acos(sqrt(zedminus));
  thetamax = M_PI - thetamin;

  const Real r2pa2 = r*r + a*a;

  gamma = E*(r2pa2*r2pa2/Delta - a*a) + a*Lz*(1. - r2pa2/Delta);
  delta = a*E*(r2pa2/Delta - 1.) - a*a*Lz/Delta;

  phi_pi2 = phi_0(0.5*M_PI);
  t_pi2 = t_0(0.5*M_PI);

  T_th = Period_Theta();
  Phi_range = Phi();

  Om_th = Omega_Theta();
  Om_phi = Omega_Phi();
  wmk = omegamk();
  pmk = wmk - m*a/(2.*Kerr::rplus(a));
}

Real CKG::Efunc()
{
  const Real rm1 = r - 1.;
  const Real numer = a*a*Lz*Lz*rm1 + r*Delta*Delta;

  //
  // This guy's a bit of a mess ...
  //
  Real denom = r*r*r*r*(r-3.);
  denom += a*a*r*(Lz*Lz + 2.*r*rm1);
  denom += a*a*a*a*(1.+r);
  denom  = Delta*sqrt(r*denom);
  denom += a*Lz*(r*r-a*a);

  return(numer/denom);
}

Real CKG::Lzfunc()
{
  const Real rm1 = r - 1.;
  const Real tmp1 = r*(1. + r*(E*E - 1.));
  const Real numer = E*(r*r - a*a) - Delta*sqrt(tmp1);
  Real denom = a*rm1;

  return(numer/denom);
}

Real CKG::Qfunc()
{
  const Real Q1 = SQR(E*(a*a + r*r) - a*Lz)/Delta;
  const Real Q2 = r*r + a*E*(a*E - 2.*Lz) + Lz*Lz;
  
  if (fabs(Q1 - Q2) < EPSILON) return(0.);
  else return(Q1 - Q2);
}

Real CKG::Omega_Theta()
{
  return(2.*M_PI/Period_Theta());
}

Real CKG::Phi()
{
  return(4.*phi_pi2);
}

Real CKG::Omega_Phi()
{
  return(Phi()/Period_Theta());
}

Real CKG::Period_Theta()
{
  return(4.*t_pi2);
}

Real CKG::dtdchi(const Real chi)
{
  const Real cchi = cos(chi);
  const Real z = zedminus*cchi*cchi;
  const Real tmp1 = gamma + a*a*E*z;
  const Real tmp2 = sqrt(betazedplus - beta_wilkins*z);

  return(tmp1/tmp2);
}

Real CKG::t(const Real chi)
{
  if (chi < 0.5*M_PI)
    return(t_pi2 - t_0(0.5*M_PI - chi));
  else
    return(t_pi2 + t_0(chi - 0.5*M_PI));
}

Real CKG::phi(const Real chi)
{
  if (chi < 0.5*M_PI)
    return(phi_pi2 - phi_0(0.5*M_PI - chi));
  else
    return(phi_pi2 + phi_0(chi - 0.5*M_PI));
}

Real CKG::t_0(const Real chi)
{
  const Real tmp1 = gamma/sqrt(betazedplus)*ellf(chi, k_zmzp);
  const Real tmp2 = (E/((1. - E)*(1. + E)))*sqrt(betazedplus)*
    (ellf(chi, k_zmzp) - elle(chi, k_zmzp));
  return(tmp1 + tmp2);
}

Real CKG::phi_0(const Real chi)
{
  //
  // Use analytic evaluation of Lz*ellpi(chi, -zedminus, k_zmzp) for
  // cosiota very close to zero.
  //
  Real term1;
  if (fabs(cosiota) < 2.e-5) {
    const Real y = beta_wilkins/Q;
    const Real z = sqrt(Q - beta_wilkins);
    const Real sign = (Lz >= 0. ? 1. : -1.);
    if (chi == 0.5*M_PI) {
      term1 = 0.5*M_PI*z*sign*(1. + 0.5*y*(1 + 0.75*y*(1. + 5.*y/6.)));
    } else {
      const Real x = Lz/sqrt(Q);
      //      
      const Real term1_0 = z*atan(Lz*tan(chi)/z);
      const Real term1_1 = 0.5*z*(atan(z/Lz) - atan(z*cos(chi)/Lz));
      const Real term1_2 = 0.375*(z*sign*atan(fabs(Lz)*tan(chi)/z) - Lz*chi);
      const Real term1_3 = 0.3125*z*(atan(z/Lz) - atan(z*cos(chi)/Lz))
	+ Lz*(cos(chi) - 1.);
      //
      term1 = term1_0 + y*(term1_1 + y*(term1_2 + y*term1_3));
    }
  } else {
    term1 = Lz*ellpi(chi, -zedminus, k_zmzp);
  }
  const Real term2 = delta*ellf(chi, k_zmzp);

  return((term1 + term2)/sqrt(betazedplus));
}

Real CKG::omegamk()
{
  return(m*Om_phi + k*Om_th);
}
