//---------------------------------------------------------------------------
//
// $Id: CEKR.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods of class CEKR --- the class containing objects relevant
// to radiation for circular, equatorial orbits.
//
// Scott Hughes, 11 January 1999.
//
#include <math.h>
#include "Globals.h"
#include "CEKR.h"
#include "SWSH.h"
#include "SNT.h"

#define SW -2

// Constructor.
CEKR::CEKR(SWSH *swsh_in, SNT *snt_in, CEKG *cekg_in) :
  swsh(swsh_in), snt(snt_in), cekg(cekg_in)
{
  //
  // Point input objects to local objects.
  //
  swsh = swsh_in;
  snt = snt_in;
  cekg = cekg_in;
  //
  // Local storage of orbital characteristics.
  //
  r = cekg->r;
  a = cekg->a;
  E = cekg->E;
  Lz = cekg->Lz;
  //
  // Local storage of mode characteristics.
  //
  l = swsh->Gl;
  m = cekg->m;
  wm = cekg->wm;
  pm = cekg->pm;
  //
  // Local storage of useful quantities that appear all over the
  // place.
  //
  Delta = Kerr::Delta(r, a);
  dr_Delta = Kerr::dr_Delta(r);
  ddr_Delta = Kerr::ddr_Delta();
  r_plus = Kerr::rplus(a);
  SN_Kay = snt->K(r);
  dr_SN_Kay = snt->dr_K(r);
  //
  // Compute amplitudes ZedH, ZedInf.
  //
  ZedH = Zed(snt->TeukRH, snt->dr_TeukRH, snt->ddr_TeukRH);
  ZedInf = Zed(snt->TeukRInf, snt->dr_TeukRInf, snt->ddr_TeukRInf);
}

// The complex amplitudes ZH, ZInf that are used to
// construct Edot, Lzdot, and the waveform.  To get
// ZH, pass in RH and its derivatives; likewise for
// ZInf.
Complex CEKR::Zed(const Complex R,
		  const Complex dr_R,
		  const Complex ddr_R)
{
  const Complex term1 = R*(Ann0() + Anmbar0() + Ambarmbar0());
  const Complex term2 = -dr_R*(Anmbar1() + Ambarmbar1());
  const Complex term3 = ddr_R*Ambarmbar2();
  const Complex denom = II*wm*snt->Bin;

  return(M_PI*(term1 + term2 + term3)/denom);
}

// The Newman-Penrose quantity rho (with sign as in
// Teukolsky's papers).
Complex CEKR::rho_func()
{
  const Complex tmp = Complex(r, 0);

  return(-1./tmp);
}

// The C_{ab} functions are projections of the orbiting
// particle's stress-energy tensor onto the Newman-Penrose
// tetrad (modulo some factors).  Here they are, plus
// their derivatives.
Complex CEKR::Cnn()
{
  const Real tmp = E*(r*r + a*a) - a*Lz;
  const Real sig = Kerr::Sigma(r, a, 0);

  return(Complex(0.25*tmp*tmp/(sig*sig*gkg.TFunc(r, a, 0., E, Lz, 0.)), 0.));
}

Complex CEKR::Cnmbar()
{
  const Real st = 1.;
  const Complex rho = rho_func();
  const Real sig = Kerr::Sigma(r, a, 0.);

  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  const Complex tmp2 = Complex(sqrt(gkg.ThetaFunc(r, a, 0., E, Lz, 0.)),
			       st*a*E - Lz/st);
  const Complex tmp3 = tmp1*tmp2;

  return(rho*tmp3/(2.*sqrt(2.)*sig*gkg.TFunc(r, a, 0., E, Lz, 0.)));
}

Complex CEKR::Cmbarmbar()
{
  const Real st = 1.;
  const Complex rho = rho_func();

  const Complex tmp = Complex(sqrt(gkg.ThetaFunc(r, a, 0., E, Lz, 0.)),
			      a*E*st - Lz/st);

  return(rho*rho*tmp*tmp/(2.*gkg.TFunc(r, a, 0., E, Lz, 0.)));
}

Complex CEKR::Ann0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;

  const Complex tmp1 = 2.*II*a*rho*st*swsh->l2dagspheroid(0.);
  const Real tmp2 = swsh->l1dagl2dagspheroid(0.);

  const Complex pref = -2.*Cnn()/(Delta*Delta*rho*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CEKR::Anmbar0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;
  const Real S = swsh->spheroid(0.);
  const Real L2dagS = swsh->l2dagspheroid(0.);

  const Complex tmp1 = (II*SN_Kay/Delta - rho - rhob)*L2dagS;
  const Complex tmp2 = (II*SN_Kay/Delta + rho +
			rhob)*II*a*st*S*(rho - rhob);

  const Complex pref = -2.*sqrt(2.)*Cnmbar()/(Delta*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CEKR::Anmbar1()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real st = 1.;
  const Real S = swsh->spheroid(0.);
  const Real L2dagS = swsh->l2dagspheroid(0.);

  const Complex tmp = L2dagS + II*a*st*(rho - rhob)*S;

  const Complex pref = -2.*sqrt(2.)*Cnmbar()/(Delta*rho3);

  return(pref*tmp);
}
    
Complex CEKR::Ambarmbar0()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  const Complex tmp1 = (SN_Kay/Delta)*(SN_Kay/Delta);
  const Complex tmp2 = 2.*rho*II*SN_Kay/Delta;
  const Complex tmp3 = II*(dr_SN_Kay - SN_Kay*dr_Delta/Delta)/Delta;

  const Complex pref = S*Cmbarmbar()*rhob/rho3;

  return(pref*(tmp1 + tmp2 + tmp3));
}
    
Complex CEKR::Ambarmbar1()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  const Complex tmp = -II*SN_Kay/Delta;

  const Complex pref = 2.*S*rhob*Cmbarmbar()/rho3;

  return(pref*(rho + tmp));
}

Complex CEKR::Ambarmbar2()
{
  const Complex rho = rho_func();
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real S = swsh->spheroid(0.);

  return(-Cmbarmbar()*rhob*S/rho3);
}

//
// Utility functions for Eric's z term.
//
Real CEKR::zerob()
{
  const Real f = Delta/(r*r);

  return(sqrt((l-1)*l*(l+1)*(l+2))*swsh->zeroY(l, m, 0.)*E/(2.*f));
}

Real CEKR::neg1b()
{
  return(sqrt((l-1)*(l+2))*swsh->neg1Y(l, m, 0.)*Lz/r);
}

Real CEKR::neg2b()
{
  return(swsh->neg2Y(l, m, 0.)*Lz*wm/((Real)m));
}

Complex CEKR::EricZed(const Complex R,
		      const Complex dr_R,
		      const Complex ddr_R)
{
  const Real zzb = zerob();
  const Real n1b = neg1b();
  const Real n2b = neg2b();
  const Real f = Delta/(r*r);
  const Complex iwrof = II*wm*r/f;

  const Complex tmp1 = zzb + 2.*II*n1b*(1. + iwrof/2.) -
    n2b*iwrof*(1. - 1./r + II*wm*r/2.)/f;

  const Complex tmp2 = -r*(II*n1b - n2b*(1. + iwrof));

  const Complex tmp3 = -n2b*r*r/2.;

  const Complex term1 = R*tmp1;
  const Complex term2 = dr_R*tmp2;
  const Complex term3 = ddr_R*tmp3;
  const Complex denom = 2.*II*wm*r*r*snt->Bin;

  return(2.0*M_PI*(term1 + term2 + term3)/denom);
}
