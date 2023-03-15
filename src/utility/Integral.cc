//---------------------------------------------------------------------------
//
// $Id: Integral.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
#include <math.h>
#include "Globals.h"
#include "Integral.h"
#include "NRUtil.h"

#define JMAX 16
#define JMAXP (JMAX+1)
#define K 5
#define TINY 1.0e-25

Complex Integral::qromb(const Real a, const Real b)
{
  int j;
  Real h[JMAXP+1];
  Complex ss, dss, s[JMAXP+1];

  h[1] = 1.0;
  for (j = 1; j <= JMAX; j++) {
    s[j] = trapzd(a, b, j);
    if (j >= K) {
      ratint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
      if (fabs(dss.imag()) <= EPSILON*fabs(ss.imag()) &&
	  fabs(dss.real()) <= EPSILON*fabs(ss.real()) &&
	  j > K + 1) return ss;
    }
    s[j+1] = s[j];
    h[j+1] = 0.25*h[j];
  }
  //
  // Value of the integral is probably zero and it's having
  // trouble converging.
  //
  cerr << "Too many steps in routine qromb..." << endl;
  cerr << "...returning ZERO" << endl;
  return 0.;
}

Complex Integral::trapzd(const Real a, const Real b, const int n)
{
  Real x, tnm, del;
  Complex sum;
  static Complex s;
  int it, j;
  
  if (n == 1) {
    return (s = 0.5*(b - a)*(FUNC(a) + FUNC(b)));
  } else {
    for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
    tnm = it;
    del = (b - a)/tnm;
    x = a + 0.5*del;
    for (sum = 0.0, j = 1; j <= it; j++, x+=del) sum += FUNC(x);
    s = 0.5*(s + (b - a)*sum/tnm);
    return s;
  }
}

void ratint(const Real xa[], const Complex ya[], const int n,
	    const Real x, Complex *y, Complex *dy)
{
  int m, i, ns=1;
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
#undef TINY
#undef K
#undef JMAXP
#undef JMAX
