//---------------------------------------------------------------------------
//
// $Id: SNTOdeint.cc,v 1.3 2004/05/29 00:55:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <math.h>
#include "Globals.h"
#include "SNT.h"
#include "NRUtil.h"

#define MAXSTP 10000000
#define TINY 1.0e-30

void SNT::odeint(Real ystart[], const int nvar, const Real x1, const Real x2,
		 const Real eps, const Real h1, const Real hmin)
{
  int nstp, i;
  Real xsav, x, hnext, hdid, h;
  Real *yscal, yscal_max = 1;
  Real *y, *dydx;

  yscal = Realvector(1, nvar);
  y = Realvector(1, nvar);
  dydx = Realvector(1, nvar);
  x = x1;
  h = SIGN(h1, x2 - x1);
  for (i = 1; i <= nvar; i++) y[i] = ystart[i];
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    SNT_derivs(x, y, dydx);
    //
    // Note the value of yscale!  This is because our amplitude and
    // derivative are typically both near zero for quite a while.  This
    // scaling is more robust.
    //
    for (i = 1; i <= nvar; i++) {
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h);
      if (0.1*yscal[i] > yscal_max) yscal_max = 0.1*yscal[i];
      yscal[i] += yscal_max;
    }
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    bsstep(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
    if ((x - x2)*(x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return;
    }
    if (fabs(hnext) <= hmin)
      Die("Step size too small in odeint");
    h = hnext;
  }
  Die("Too many steps in routine odeint");
}

void SNT::odeint(Complex ystart[], const int nvar, const Real x1,
		 const Real x2, const Real eps, const Real h1,
		 const Real hmin)
{
  int nstp, i;
  Real xsav, x, hnext, hdid, h;
  Real *yscal;
  Complex *y, *dydx;
  
  yscal = Realvector(1, nvar);
  y = Complexvector(1, nvar);
  dydx = Complexvector(1, nvar);
  x = x1;
  h = SIGN(h1, x2 - x1);
  for (i = 1; i <= nvar; i++) y[i] = ystart[i];
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    SNT_derivs(x, y, dydx);
    for (i = 1; i <= nvar; i++)
      yscal[i] = FMAX(fabs( y[i].real() ) + fabs( dydx[i].real()*h ) + TINY,
		     fabs( y[i].imag() ) + fabs( dydx[i].imag()*h ) + TINY);
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    bsstep(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
    if ((x - x2)*(x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      free_Complexvector(dydx, 1, nvar);
      free_Complexvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return;
    }
    if (fabs(hnext) <= hmin) Die("Step size too small in odeint");
    h = hnext;
  }
  Die("Too many steps in routine odeint");
}

#undef MAXSTP
#undef TINY
