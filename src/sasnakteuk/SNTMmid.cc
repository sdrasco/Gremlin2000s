//---------------------------------------------------------------------------
//
// $Id: SNTMmid.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include "Globals.h"
#include "SNT.h"
#include "NRUtil.h"

void SNT::mmid(Real y[], Real dydx[], int nvar, Real xs, Real htot,
	       int nstep, Real yout[])
{
  int n, i;
  Real swap, *ym, *yn;
  Real x, h2, h;
  
  ym = Realvector(1, nvar);
  yn = Realvector(1, nvar);
  h = htot/nstep;
  for (i = 1; i <= nvar; i++) {
    ym[i] = y[i];
    yn[i] = y[i] + h*dydx[i];
  }
  x = xs + h;
  SNT_derivs(x, yn, yout);

  h2 = 2.0*h;
  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    SNT_derivs(x, yn, yout);
  }
  for (i = 1; i<= nvar; i++)
    yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
  free_Realvector(yn, 1, nvar);
  free_Realvector(ym, 1, nvar);
}

void SNT::mmid(Complex y[], Complex dydx[], int nvar, Real xs,
	       Real htot, int nstep, Complex yout[])
{
  int n, i;
  Complex swap, *ym, *yn;
  Real x, h2, h;
  
  ym = Complexvector(1, nvar);
  yn = Complexvector(1, nvar);
  h = htot/nstep;
  for (i = 1; i <= nvar; i++) {
    ym[i] = y[i];
    yn[i] = y[i] + h*dydx[i];
  }
  x = xs + h;
  SNT_derivs(x, yn, yout);

  h2 = 2.0*h;
  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    SNT_derivs(x, yn, yout);
  }
  for (i = 1; i<= nvar; i++)
    yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
  free_Complexvector(yn, 1, nvar);
  free_Complexvector(ym, 1, nvar);
}
