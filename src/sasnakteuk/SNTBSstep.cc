//---------------------------------------------------------------------------
//
// $Id: SNTBSstep.cc,v 1.2 2004/09/28 16:05:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <math.h>
#include "Globals.h"
#include "SNT.h"
#include "NRUtil.h"
#include "NRSNT.h"

#define KMAXX 8
#define IMAXX (KMAXX + 1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

void SNT::bsstep(Real y[], Real dydx[], int nv, Real *xx, Real htry,
		 Real eps, Real yscal[], Real *hdid, Real *hnext)
{
  int i, iq, k, kk, km;
  static int first = 1, kmax, kopt;
  static Real epsold = -1.0, xnew;
  Real **d;
  Real *x;
  Real eps1, errmax, fact, h, red, scale, work, wrkmin, xest;
  Real *err;
  Real *yerr, *ysav, *yseq;
  static Real a[IMAXX + 1];
  static Real alf[KMAXX + 1][KMAXX + 1];
  static int nseq[IMAXX + 1] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
  int reduct, exitflag = 0;

  d = Realmatrix(1, KMAXX, 1, KMAXX);
  err = Realvector(1, KMAXX);
  x = Realvector(1, KMAXX);
  yerr = Realvector(1, nv);
  ysav = Realvector(1, nv);
  yseq = Realvector(1, nv);
  if (eps != epsold) {
    *hnext = xnew = -1.0e29;
    eps1 = SAFE1*eps;
    a[1] = nseq[1] + 1;
    for (k = 1; k <= KMAXX; k++) a[k + 1] = a[k] + nseq[k + 1];
    for (iq = 2; iq <= KMAXX; iq++) {
      for (k = 1; k < iq; k++)
	alf[k][iq] = pow(eps1, (a[k + 1] - a[iq + 1])/
			 ((a[iq + 1] - a[1] + 1.0)*(2*k + 1)));
    }
    epsold = eps;
    for (kopt = 2; kopt < KMAXX; kopt++)
      if (a[kopt + 1] > a[kopt]*alf[kopt - 1][kopt]) break;
    kmax = kopt;
  }
  h = htry;
  for (i = 1; i <= nv; i++) ysav[i] = y[i];
  if (*xx != xnew || h != (*hnext)) {
    first = 1;
    kopt = kmax;
  }
  reduct = 0;
  for (;;) {
    for (k = 1; k <= kmax; k++) {
      xnew = (*xx) + h;
      if (xnew == (*xx)) Die("step size underflow in bsstep");
      mmid(ysav, dydx, nv, *xx, h, nseq[k], yseq);
      xest = SQR(h/nseq[k]);
      rzextr(k, xest, yseq, y, yerr, d, x, nv);
      if (k != 1) {
	errmax = TINY;
	for (i = 1; i <= nv; i++)
	  errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	km = k - 1;
	err[km] = pow(errmax/SAFE1, 1.0/(2*km + 1));
      }
      if (k != 1 && (k >= kopt - 1 || first)) {
	if (errmax < 1.0) {
	  exitflag = 1;
	  break;
	}
	if (k == kmax || k == kopt + 1) {
	  red = SAFE2/err[km];
	  break;
	}
	else if (k == kopt && alf[kopt - 1][kopt] < err[km]) {
	  red = 1.0/err[km];
	  break;
	}
	else if (kopt == kmax && alf[km][kmax - 1] < err[km]) {
	  red = alf[km][kmax - 1]*SAFE2/err[km];
	  break;
	}
	else if (alf[km][kopt] < err[km]) {
	  red = alf[km][kopt - 1]/err[km];
	  break;
	}
      }
    }
    if (exitflag) break;
    red = FMIN(red, REDMIN);
    red = FMAX(red, REDMAX);
    h *= red;
    reduct = 1;
  }
  *xx = xnew;
  *hdid = h;
  first = 0;
  wrkmin = 1.0e35;
  for (kk = 1; kk <= km; kk++) {
    fact = FMAX(err[kk], SCALMX);
    work = fact*a[kk + 1];
    if (work < wrkmin) {
      scale = fact;
      wrkmin = work;
      kopt = kk + 1;
    }
  }
  *hnext = h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact = FMAX(scale/alf[kopt - 1][kopt], SCALMX);
    if (a[kopt + 1]*fact <= wrkmin) {
      *hnext = h/fact;
      kopt++;
    }
  }
  free_Realvector(yseq, 1, nv);
  free_Realvector(ysav, 1, nv);
  free_Realvector(yerr, 1, nv);
  free_Realvector(x, 1, KMAXX);
  free_Realvector(err, 1, KMAXX);
  free_Realmatrix(d, 1, KMAXX, 1, KMAXX);
}

void SNT::bsstep(Complex y[], Complex dydx[], int nv, Real *xx,
		 Real htry, Real eps, Real yscal[], Real *hdid,
		 Real *hnext)
{
  int i, iq, k, kk, km;
  static int first = 1, kmax, kopt;
  static Real epsold = -1.0, xnew;
  Complex **d;
  Real *x;
  Real eps1, errmax, fact, h, red, scale, work, wrkmin, xest;
  Real *err;
  Complex *yerr, *ysav, *yseq;
  static Real a[IMAXX + 1];
  static Real alf[KMAXX + 1][KMAXX + 1];
  static int nseq[IMAXX + 1] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
  int reduct, exitflag = 0;

  d = Complexmatrix(1, KMAXX, 1, KMAXX);
  err = Realvector(1, KMAXX);
  x = Realvector(1, KMAXX);
  yerr = Complexvector(1, nv);
  ysav = Complexvector(1, nv);
  yseq = Complexvector(1, nv);
  if (eps != epsold) {
    *hnext = xnew = -1.0e29;
    eps1 = SAFE1*eps;
    a[1] = nseq[1] + 1;
    for (k = 1; k <= KMAXX; k++) a[k + 1] = a[k] + nseq[k + 1];
    for (iq = 2; iq <= KMAXX; iq++) {
      for (k = 1; k < iq; k++)
	alf[k][iq] = pow(eps1, (a[k + 1] - a[iq + 1])/
		       ((a[iq + 1] - a[1] + 1.0)*(2*k + 1)));
    }
    epsold = eps;
    for (kopt = 2; kopt < KMAXX; kopt++)
      if (a[kopt + 1] > a[kopt]*alf[kopt - 1][kopt]) break;
    kmax = kopt;
  }
  h = htry;
  for (i = 1; i <= nv; i++) ysav[i] = y[i];
  if (*xx != xnew || h != (*hnext)) {
    first = 1;
    kopt = kmax;
  }
  reduct = 0;
  for (;;) {
    for (k = 1; k <= kmax; k++) {
      xnew = (*xx) + h;
      if (xnew == (*xx)) Die("step size underflow in SNT::bsstep(Complex y[] ... )");
      mmid(ysav, dydx, nv, *xx, h, nseq[k], yseq);
      xest = SQR(h/nseq[k]);
      rzextr(k, xest, yseq, y, yerr, d, x, nv);
      if (k != 1) {
	errmax = TINY;
	for (i = 1; i <= nv; i++)
	  errmax = FMAX(errmax, abs(yerr[i]/yscal[i]));
	errmax /= eps;
	km = k - 1;
	err[km] = pow(errmax/SAFE1, 1.0/(2*km + 1));
      }
      if (k != 1 && (k >= kopt - 1 || first)) {
	if (errmax < 1.0) {
	  exitflag = 1;
	  break;
	}
	if (k == kmax || k == kopt + 1) {
	  red = SAFE2/err[km];
	  break;
	}
	else if (k == kopt && alf[kopt - 1][kopt] < err[km]) {
	  red = 1.0/err[km];
	  break;
	}
	else if (kopt == kmax && alf[km][kmax - 1] < err[km]) {
	  red = alf[km][kmax - 1]*SAFE2/err[km];
	  break;
	}
	else if (alf[km][kopt] < err[km]) {
	  red = alf[km][kopt - 1]/err[km];
	  break;
	}
      }
    }
    if (exitflag) break;
    red = FMIN(red, REDMIN);
    red = FMAX(red, REDMAX);
    h *= red;
    reduct = 1;
  }
  *xx = xnew;
  *hdid = h;
  first = 0;
  wrkmin = 1.0e35;
  for (kk = 1; kk <= km; kk++) {
    fact = FMAX(err[kk], SCALMX);
    work = fact*a[kk + 1];
    if (work < wrkmin) {
      scale = fact;
      wrkmin = work;
      kopt = kk + 1;
    }
  }
  *hnext = h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact = FMAX(scale/alf[kopt - 1][kopt], SCALMX);
    if (a[kopt + 1]*fact <= wrkmin) {
      *hnext = h/fact;
      kopt++;
    }
  }
  free_Complexvector(yseq, 1, nv);
  free_Complexvector(ysav, 1, nv);
  free_Complexvector(yerr, 1, nv);
  free_Realvector(x, 1, KMAXX);
  free_Realvector(err, 1, KMAXX);
  free_Complexmatrix(d, 1, KMAXX, 1, KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
