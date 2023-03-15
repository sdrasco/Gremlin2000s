//---------------------------------------------------------------------------
//
// $Id: SNTTeuk.cc,v 1.4 2004/06/01 14:15:17 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Code to transform the Sasaki-Nakamura solution to the Teukolsky solution.
//
// Scott Hughes, 9 January 1999
//
#include <math.h>
#include "Globals.h"
#include "SNT.h"

Complex SNT::CalcR(const Complex X, const Complex dr_X, const Real rad)
{
  const Real Delta = Kerr::Delta(rad, a);
  const Real dr_Delta = Kerr::dr_Delta(rad);
  const Real r2pa2 = rad*rad + a*a;
  const Real sr2pa2 = sqrt(r2pa2);
  
  const Complex chi = X*Delta/sr2pa2;
  const Complex dr_chi = dr_X*Delta/sr2pa2 +
    X*(dr_Delta - rad*Delta/r2pa2)/sr2pa2;
  
  return((chi*(alpha(rad) + dr_beta(rad)/Delta) -
	  dr_chi*beta(rad)/Delta)/eta(rad));
}

Complex SNT::Calcdr_R(const Complex X, const Complex dr_X,
		      const Complex ddr_X, const Real rad)
{
  const Real Delta = Kerr::Delta(rad, a);
  const Real dr_Delta = Kerr::dr_Delta(rad);
  const Real ddr_Delta = Kerr::ddr_Delta();
  const Real r2pa2 = rad*rad + a*a;
  const Real sr2pa2 = sqrt(r2pa2);
  
  const Complex chi = X*Delta/sr2pa2;
  const Complex dr_chi = dr_X*Delta/sr2pa2 +
    X*(dr_Delta - rad*Delta/r2pa2)/sr2pa2;
  const Complex ddr_chi = ddr_X*Delta/sr2pa2 +
    2.*dr_X*(dr_Delta - rad*Delta/r2pa2)/sr2pa2 +
    X*(ddr_Delta - (2.*rad*dr_Delta + Delta)/r2pa2 +
       3.*rad*rad*Delta/(r2pa2*r2pa2))/sr2pa2;

  const Complex alp = alpha(rad);
  const Complex dr_alp = dr_alpha(rad);
  const Complex bet = beta(rad);
  const Complex dr_bet = dr_beta(rad);
  const Complex ddr_bet = ddr_beta(rad);
  const Complex eeta = eta(rad);

  const Complex ans1 = -dr_eta(rad)*(chi*(alp + dr_bet/Delta) -
				     dr_chi*bet/Delta)/(eeta*eeta);

  const Complex ans2 = ((dr_alp + ddr_bet/Delta -
			 dr_bet*dr_Delta/(Delta*Delta))*chi +
			(alp + bet*dr_Delta/(Delta*Delta))*dr_chi -
			(bet/Delta)*ddr_chi)/eeta;
  return(ans1 + ans2);
}

Complex SNT::Calcddr_R(Complex R, Complex dr_R, const Real rad)
{
  const Real Delta = Kerr::Delta(rad, a);
  const Real dr_Delta = Kerr::dr_Delta(rad);

  return((dr_Delta/Delta)*dr_R + (V(rad)/Delta)*R);
}
