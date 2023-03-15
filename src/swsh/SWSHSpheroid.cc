//---------------------------------------------------------------------------
//
// $Id: SWSHSpheroid.cc,v 1.3 2004/09/17 19:33:06 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// The functions in this file are used to compute the spin weighted
// spheroidal harmonics.  See my notes for mathematical details.
// These functions are used by the class SWSH.
//
// Scott Hughes, 24 July 1998
//

#include <math.h>
#include "Globals.h"
#include "SWSH.h"
#include "NRUtil.h"
#include "NRSWSH.h"

SWSH::SWSH(const int ll, const int ss, const int mm,
	   const Real spintimesfreq) :
  Gl(ll), s(ss), Gm(mm), aw(spintimesfreq)
{
  lmin = IMAX(abs(s), abs(Gm));

  expand(&E, b, &N);
  lambda = E - 2.*Gm*aw + aw*aw - 2.;

  cheb_spheroid();
  cheb_l2dagspheroid();
  cheb_l1dagl2dagspheroid();
}

//
// Wrapper function for load_M and findeig.  This function starts out
// keeping a small number of terms in b, and then allows b to grow until
// findeig says that the b vector is OK.
//
void SWSH::expand(Real *E, Real *b, int *n)
{
  Real **M, norm = 0.0, maxcof;
  int okflag = 0, i, K, maxind;

  if (Gl < abs(s) || Gl < abs(Gm))
    Die("Illegal values of m or s in expand");
  K = Gl - lmin + 7;
  while(!okflag) {
    M = Realmatrix(1, K, 1, K);
    load_M(M, K);
    okflag = findeig(M, K, E, b);
    if (!okflag)
      free_Realmatrix(M, 1, K, 1, K);
    K += 2;
  }
  K -= 2;
  maxcof = 0.0;
  for (i = 1; i <= K; i++) {
    norm += b[i]*b[i];
    if (fabs(b[i]) > maxcof) {
      maxcof = fabs(b[i]);
      maxind = i;
    }
  }
  norm = sqrt(norm);
  if (b[maxind] < 0.0) norm *= -1.0;
  for (i = 1; i <= K; i++)
    b[i] /= norm;
  *n = K;
  free_Realmatrix(M, 1, K, 1, K);
}

void SWSH::load_M(Real **M, const int K)
{
  Clebsch cg;
  int i, j, p, q;

  for (i = 1; i <= K; i++) {
    p = i - 1 + lmin;
    for (j = 1; j <= K; j++) {
      q = j - 1 + lmin;
      if (i == j) {
	M[i][j] = aw*aw*cg.xsqrbrac(s, p, q, Gm) - 
	  2.0*aw*((Real)s)*cg.xbrac(s, p, q, Gm)
	   - ((Real)(p*(p + 1)));
      }
      else if (i == j - 1 || i == j + 1)
	M[i][j] = aw*aw*cg.xsqrbrac(s, p, q, Gm) - 
	  2.0*aw*((Real)s)*cg.xbrac(s, p, q, Gm);
      else if (i == j - 2 || i == j + 2)
	M[i][j] = aw*aw*cg.xsqrbrac(s, p, q, Gm);
      else M[i][j] = 0.0;
    }
  }
}

//
// This function finds the eigenvalues and eigenvectors of M;
// see Numerical Recipes for a description of how the routines
// tred2() (which performs a Householder transformation to put M
// in tridiagonal form) and tqli() (which then computes the
// eigenquantities) work.
//
// Housekeeping note: it's unclear a priori how many terms need
// to be kept in the vector b.  Note that this function is integer:
// it returns 1 if the last few entries are less than 10^(-15), and
// zero otherwise.  This allows the calling function to determine
// if it is keeping enough terms.
//

int SWSH::findeig(Real **M, const int K, Real *E, Real b[])
{
  Real *e, *d;
  int i;
  long unsigned int *ind;
  
  e = Realvector(1, K);
  d = Realvector(1, K);
  ind = (long unsigned int*)ivector(1, K);
  for (i = 1; i <= K; i++) ind[i] = i;
  tred2(M, K, d, e);
  tqli(d, e, K, M);
  for (i = 1; i <= K; i++) d[i] = -d[i];
  if (s == 0) indexx(K, d, ind);
  if (fabs(M[K][Gl + 1 - lmin]) < 1.e-15
     && fabs(M[K - 1][Gl + 1 - lmin]) < 1.e-15) {
    *E = d[ind[Gl - lmin + 1]];
    for (i = 1; i <= K; i++)
      b[i] = M[i][ind[Gl + 1 - lmin]];
    free_Realvector(e, 1, K);
    free_Realvector(d, 1, K);
    free_ivector((int *)ind, 1, K);
    return(1);
  } else {
    free_Realvector(e, 1, K);
    free_Realvector(d, 1, K);
    free_ivector((int *)ind, 1, K);
    return(0);
  }
}

Real SWSH::error(const Real E, const Real b[], const int n)
{
  Real **M, *berr;
  Real err = 0.0;
  int i, j;

  M = Realmatrix(1, n, 1, n);
  load_M(M, n);
  berr = Realvector(1, n);
  for (i = 1; i <= n; i++) {
    berr[i] = 0.0;
    for (j = 1; j <= n; j++)
      berr[i] += M[i][j]*b[j];
    berr[i] += E*b[i];
  }
  for (i = 1; i <= n; i++)
    err += berr[i]*berr[i];
  free_Realvector(berr, 1, n);
  free_Realmatrix(M, 1, n, 1, n);
  return(sqrt(err)/E);
}

//
// Returns spheroid.
//
Real SWSH::spheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  y2 = 2.0*(y = x);
  for (j = N_spheroid; j >= 1; j--) {
    sv = d;
    d = y2*d - dd + spheroid_cofs[j];
    dd = sv;
  }
  if (Gm%2)
    ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*spheroid_cofs[0]);
  else
    ans = y*d - dd + 0.5*spheroid_cofs[0];
  return ans;
}

//
// Returns spheroid.
//
Real SWSH::l2dagspheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  y2 = 2.0*(y = x);
  for (j = N_l2dagspheroid; j >= 1; j--) {
    sv = d;
    d = y2*d - dd + l2dagspheroid_cofs[j];
    dd = sv;
  }
  if (Gm%2)
    ans = y*d - dd + 0.5*l2dagspheroid_cofs[0];
  else
    ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*l2dagspheroid_cofs[0]);
  return ans;
}

//
// Returns spheroid.
//
Real SWSH::l1dagl2dagspheroid(const Real x)
{
  Real d = 0.0, dd = 0.0, sv, y, y2;
  Real ans;
  int j;

  y2 = 2.0*(y = x);
  for (j = N_l1dagl2dagspheroid; j >= 1; j--) {
    sv = d;
    d = y2*d - dd + l1dagl2dagspheroid_cofs[j];
    dd = sv;
  }
  if (Gm%2)
    ans = sqrt((1. - x)*(1. + x))*(y*d - dd + 0.5*l1dagl2dagspheroid_cofs[0]);
  else
    ans = y*d - dd + 0.5*l1dagl2dagspheroid_cofs[0];
  return ans;
}

void SWSH::cheb_spheroid()
{
  int k;
  Real y, fac;

  fac = 2.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      spheroid_cofs[k] = fac*raw_spheroid(y)/sqrt((1. - y)*(1. + y));
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      spheroid_cofs[k] = fac*raw_spheroid(y);
    }
  }
  cosft2(spheroid_cofs - 1, CHEBCOFS);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(spheroid_cofs[k]) > 5.e-15) N_spheroid = k;
}

void SWSH::cheb_l2dagspheroid()
{
  int k;
  Real y, fac;

  fac = 2.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      l2dagspheroid_cofs[k] = fac*raw_l2dagspheroid(y);
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      l2dagspheroid_cofs[k] = fac*raw_l2dagspheroid(y)/sqrt((1. - y)*(1. + y));
    }
  }
  cosft2(l2dagspheroid_cofs - 1, CHEBCOFS);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(l2dagspheroid_cofs[k]) > 5.e-15) N_l2dagspheroid = k;
}

void SWSH::cheb_l1dagl2dagspheroid()
{
  int k;
  Real y, fac;

  fac = 2.0/CHEBCOFS;
  if (Gm%2) {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      l1dagl2dagspheroid_cofs[k] =
	fac*raw_l1dagl2dagspheroid(y)/sqrt((1. - y)*(1. + y));
    }
  } else {
    for (k = 0; k < CHEBCOFS; k++) {
      y = cos(M_PI*(k + 0.5)/CHEBCOFS);
      l1dagl2dagspheroid_cofs[k] = fac*raw_l1dagl2dagspheroid(y);
    }
  }
  cosft2(l1dagl2dagspheroid_cofs - 1, CHEBCOFS);
  for (k = 0; k < CHEBCOFS; k++)
    if (fabs(l1dagl2dagspheroid_cofs[k]) > 5.e-15) N_l1dagl2dagspheroid = k;
}

Real SWSH::raw_spheroid(const Real x)
{
  int k;
  Real ans = 0.0;

  if (s != 0 && s != -1 && s != -2)
    Die("This code can only handle spin weights 0, -1, -2.");

  for (k = lmin; k <= lmin + N - 1; k++) {
    if (s == 0)
      ans += b[k - lmin + 1]*zeroY(k, Gm, x);
    else if (s == -1)
      ans += b[k - lmin + 1]*neg1Y(k, Gm, x);
    else if (s == -2)
      ans += b[k - lmin + 1]*neg2Y(k, Gm, x);
  }
  return ans;
}

Real SWSH::raw_l2dagspheroid(const Real x)
{
  int k, s = -2;
  Real ans = 0.0;

  for (k = lmin; k <= lmin + N - 1; k++) {
    ans += aw*sqrt((1. - x)*(1. + x))*b[k - lmin + 1]*neg2Y(k, Gm, x);
    ans -= b[k - lmin + 1]*sqrt((k - 1)*(k + 2))*neg1Y(k, Gm, x);
  }

  return ans;
}

Real SWSH::raw_l1dagl2dagspheroid(const Real x)
{
  int k, s = -2;
  Real ans = 0.0;

  for (k = lmin; k <= lmin + N - 1; k++) {
    ans += sqrt((k - 1)*k*(k + 1)*(k + 2))*b[k - lmin + 1]*zeroY(k, Gm, x);
    ans -= aw*aw*(1. - x)*(1. + x)*b[k - lmin + 1]*neg2Y(k, Gm, x);
    //
    // These next two lines are the $2 aw sin(theta) L_2^\dag S$ term
    //
    ans += 2.0*aw*aw*(1. - x)*(1. + x)*b[k - lmin + 1]*neg2Y(k, Gm, x);
    ans -= 2.0*aw*sqrt((1. - x)*(1. + x))*b[k - lmin + 1]*
      sqrt((k - 1)*(k + 2))*neg1Y(k, Gm, x);
  }
  return ans;
}

//
// Only does forward transform, since that is all that is needed
// for the Chebyshev decomposition.
//
void SWSH::cosft2(Real y[], int n)
{
  int i;
  Real sum, sum1, y1, y2, ytemp;
  Real theta, wi = 0.0, wi1, wpi, wpr, wr = 1.0, wr1, wtemp;

  theta = 0.5*M_PI/n;
  wr1 = cos(theta);
  wi1 = sin(theta);
  wpr = -2.0*wi1*wi1;
  wpi = sin(2.0*theta);
  for (i = 1; i <= n/2; i++) {
    y1 = 0.5*(y[i] + y[n - i + 1]);
    y2 = wi1*(y[i] - y[n - i + 1]);
    y[i] = y1 + y2;
    y[n - i + 1] = y1 - y2;
    wr1 = (wtemp = wr1)*wpr - wi1*wpi + wr1;
    wi1 = wi1*wpr + wtemp*wpi + wi1;
  }
  realft(y, n, 1);
  for (i = 3; i <= n; i += 2) {
    wr = (wtemp = wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
    y1 = y[i]*wr - y[i + 1]*wi;
    y2 = y[i + 1]*wr + y[i]*wi;
    y[i] = y1;
    y[i + 1] = y2;
  }
  sum = 0.5*y[2];
  for (i = n; i >= 2; i -= 2) {
    sum1 = sum;
    sum += y[i];
    y[i] = sum1;
  }
}

void SWSH::realft(Real data[], unsigned long n, int isign)
{
  unsigned long i, i1, i2, i3, i4, np3;
  Real c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  Real wr, wi, wpr, wpi, wtemp, theta;

  theta = M_PI/(Real) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, n >> 1, 1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n + 3;
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

void SWSH::four1(Real data[], unsigned long nn, int isign)
{
  unsigned long n, mmax, m, j, istep, i;
  Real wtemp, wr, wpr, wpi, wi, theta;
  Real tempr, tempi;

  n = nn << 1;
  j = 1;
  for (i = 1; i<n; i += 2) {
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
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
	j = i + mmax;
	tempr = wr*data[j] - wi*data[j + 1];
	tempi = wr*data[j + 1] + wi*data[j];
	data[j] = data[i] - tempr;
	data[j + 1] = data[i + 1] - tempi;
	data[i] += tempr;
	data[i + 1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}
