//---------------------------------------------------------------------------
//
// $Id: SWSHCGUtil.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// The functions in this file are used to calculate various quantities
// related to Clebsch-Gordan coefficients.  They are used by the class
// Clebsch.
//
// Scott Hughes, 24 July 1998
//

#include <math.h>
#include "Globals.h"
#include "NRUtil.h"
#include "SWSH.h"
#include "NRSWSH.h"

//
// As background to these functions, let |s,l,m> represent the spherical
// (NOT spheroidal!) harmonic of spin weight s and with indices l and m.
// Further, let
//
// <s, p, m | f(\theta,\phi) | s, q, m>
//
// represent the integral
//
// $\int d\Omega {_sY^m_p}^* f(\theta,\phi) {_sY^m_q}
//

// 
// This function computes the integral
//
// <s, p, m | \cos\theta | s, q, m>
//
// See my notes on spin-weighted spheroidal harmonics for the
// formulae used to generate this.
//

Real Clebsch::xbrac(const int s, const int q, const int p, const int m)
{
  Real ans;

  ans = cgcof(q, 1, m, 0, p, m)*cgcof(q, 1, -s, 0, p, -s);
  ans *= sqrt(((Real)(2*q + 1))/((Real)(2*p + 1)));
  return(ans);
}

// 
// This function computes the integral
//
// <s,p,m | (\cos\theta)^2 | s,q,m>
//
// See my notes on spin-weighted spheroidal harmonics for the
// formulae used to generate this.
//

Real Clebsch::xsqrbrac(const int s, const int q, const int p, const int m)
{
  Real ans;

  ans = cgcof(q, 2, m, 0, p, m)*cgcof(q, 2, -s, 0, p, -s);
  ans *= (2.0/3.0)*sqrt(((Real)(2*q + 1))/((Real)(2*p + 1)));
  if (q == p) ans += (1.0/3.0);
  return(ans);
}

// 
// This function computes the integral
//
// <s,p,m | \sin\theta | s,q,mp>
//
// See Marc's and my notes for the formulae used to generate this.
//

Real Clebsch::sinthetabrac(const int s, const int p, const int q,
			   const int m, const int mp)
{
  Real ans;

  if (mp != m + 1 && mp != m - 1) {
    ans = 0.;
  } else if (mp == m + 1) {
    ans = cgcof(q, 1, mp, -1, p, m)*cgcof(q, 1, -s, 0, p, -s);
    ans *= sqrt(((Real)(4*q + 2))/((Real)(2*p + 1)));
  } else { // mp == m - 1 only possibility left ...
    ans = cgcof(q, 1, mp, 1, p, m)*cgcof(q, 1, -s, 0, p, -s);
    ans *= -sqrt(((Real)(4*q + 2))/((Real)(2*p + 1)));
  }
  return(ans);
}

//
// This function computes the Clebsch-Gordan coefficient:
// <j1, j2, m1, m2 | J, M>.  The algorithm here is from the textbook
// by Hamermesh ("Group theory and its application to physical
// problems").
//

Real Clebsch::cgcof(const int j1, const int j2, const int m1, const int m2,
		    const int J, const int M)
{
  int A, B, C, D, E;
  int lamb;
  Real rho, CJ = 0.0, term;

  if (m1 + m2 != M) return(0.0);
  if ((abs(j1 - j2) > J) || (j1 + j2) < J) return(0.0);
  for (lamb = 0; lamb <= (j1 + j2 - J); lamb++) {
    A = j1 + j2 - J - lamb;
    B = j1 - lamb - m1;
    C = J - j2 + lamb + m1;
    D = j2 - lamb + m2;
    E = J - j1 + lamb - m2;
    if (A < 0 || B < 0 || C < 0 || D < 0 || E < 0) 
      term = 0.0;
    else {
      if (lamb & 1) term = -1.0;
      else term = 1.0;
      term *= sqrt(factrl(j1 + m1)*factrl(j1 - m1)*factrl(j2 + m2)*
		   factrl(j2 - m2));
      term /= factrl(lamb)*factrl(A)*factrl(B);
      term *= sqrt(factrl(J + M)*factrl(J - M));
      term /= factrl(C)*factrl(D)*factrl(E);
    }
    CJ += term;
  }
  rho = sqrt(((Real)2*J + 1)*factrl(J - j2 + j1)*factrl(J - j1 + j2)*
	     factrl(j1 + j2 - J));
  rho /= sqrt(factrl(j1 + j2 + J + 1));
  return(rho*CJ);
}
