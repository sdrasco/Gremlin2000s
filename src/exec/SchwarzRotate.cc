//  Input:  e p theta_inc l m k n
//
//  output:  (E*(theta_inc) - E*(0) |D|^2 ) / E*(theta_inc) 
//
// By: Steve Drasco

// standard includes
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>

// Scott's includes
#include "CKG.h"
#include "SNT.h"
#include "SWSH.h"
#include "Globals.h"
#include "CKR.h"

// Steve's includes
#include "IEKG.h"
#include "IEKR.h"

Real factrl(int n)
{
  if (n < 0) exit(0);

  Real output = 1.0;
  for(int j = n; j > 0; j--) {
    output *= j;
  }
  return output;
}


Real wigner_d_arfken(const int l, const int m, const int mprime, const Real theta)
{
  Real factrl(int n);
  Real numer, denom, angles, sign, sum = 0.;

  // terms which don't depend on summation index
  numer = sqrt(factrl(l + mprime) * factrl(l - mprime) * factrl(l + m) * factrl(l - m));
  Real sto2 = sin(0.5*theta);
  Real cto2 = cos(0.5*theta);

  // compute sum boundaries
  int qfirst, qlast;
  if(mprime - m > 0) {
    qfirst = mprime - m;
  } else {
    qfirst = 0;
  }
  if(l - m < l + mprime) {
    qlast = l - m;
  } else {
    qlast = l + mprime;
  }

  // main sum
  for (int q = qfirst; q <= qlast; q++) {
    if (q%2) sign = -1;
    else sign = 1;
    denom = factrl(q) * factrl(l - m - q) * factrl(l + mprime - q) * factrl(m - mprime + q);
    angles = pow(cto2, Real(2*l + mprime - m - 2*q)) * pow(-sto2, Real(m - mprime + 2*q));
    sum += sign*numer*angles/denom;
  }
  return(sum);
}


Real wigner_d_www(const int l, const int m, const int mprime, const Real theta)
{
  Real factrl(int n);
  Real numer, denom, angles, sign, sum = 0.;

  // terms which don't depend on summation index
  if (l+mprime%2){
    sign = -1;
  } else {
    sign = 1;
  }
  numer = sign*sqrt(factrl(l + m) * factrl(l - m) * factrl(l + mprime) * factrl(l - mprime));
  const Real sto2 = sin(0.5*theta);
  const Real cto2 = cos(0.5*theta);

  // compute sum boundaries
  int kfirst, klast;
  if(mprime + m > 0) {
    kfirst = mprime + m;
  } else {
    kfirst = 0;
  }
  if(l + mprime < l + m) {
    klast = l+mprime;
  } else {
    klast = l+m;
  }

  // main sum
  for (int k = kfirst; k <= klast; k++) {
    if (k%2) sign = -1;
    else sign = 1;
    denom = factrl(k) * factrl(l + mprime - k) * factrl(l + m - k) * factrl(k - mprime - m);
    angles = pow(cto2, Real(2*k - m - mprime)) * pow(sto2, Real(2*l + mprime + m - 2*k));
    sum += sign*numer*angles/denom;
  }
  return(sum);
}

Real wigner_d_old(const int l, const int m, const int mprime, const Real theta)
{
  Real factrl(int n);
  Real numer, denom, angles, sign, sum = 0.;
  int k;

  for (k = 0; k <= l + m; k++) {
    if (k%2) sign = -1;
    else sign = 1;

    numer = sqrt(factrl(l + m)*factrl(l - m)*
                 factrl(l + mprime)*factrl(l - mprime));

    if (l - mprime - k < 0 || mprime - m + k < 0) {
      denom = 1;
      sign = 0.0;
    } else {
      denom = factrl(k)*factrl(l - mprime - k)*
        factrl(l + m - k)*factrl(mprime - m + k);
    }
    const Real sto2 = sin(0.5*theta);
    const Real cto2 = cos(0.5*theta);
    angles = pow(cto2, (Real)(2*l + m - mprime - 2*k))*
      pow(-sto2, (Real)(mprime - m + 2*k));

    sum += sign*numer*angles/denom;
  }
  return(sum);
}

void GetResidual(Real e, Real p, Real theta_inc, int l, int m, int k, int n, Real eps, Real &circular_resid, Real &generic_resid)
{

  // compute theta_
  theta_inc *= M_PI / 180;
  Real theta_;
  int prograde;
  if(theta_inc < 0.5 * M_PI) {
    prograde = 1;
    theta_ = 0.5*M_PI - theta_inc;
  } else {
    prograde = 0;
    theta_ = theta_inc - 0.5*M_PI;
  }

  // Schwarzschild limit
  Real a = 0.0;

  // equatorial orbit
  IEKG iekg(prograde,a,e,p,0.5*M_PI,Min(eps, 1e-12));
  IEKR iekr(&iekg,l,m,0,n,eps); 
  Complex ZH = iekr.ZH;
  Complex ZInf = iekr.ZInf;
  Real alpha = iekr.alpha;
  Real omega = iekg.Omega(m,0,n);
  Real EInf = ZH.real()*ZH.real() + ZH.imag()*ZH.imag();
  EInf /= 4.0*M_PI*omega*omega;
  Real EH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
  EH *= alpha;
  EH /= 4.0*M_PI*omega*omega;

  // equatorial circular orbit
  Real MinR = p/(1.0 + e);
  Real MaxR = p/(1.0 - e);
  Real AverageR = 0.5*(MinR+MaxR);
  Real r = AverageR;
  Real circ_a = 1e-4;
  CKG ckg(m, 0, iekg.CosIota, USING_cosiota, r, circ_a);
  SWSH swsh_circ(l, -2, m, circ_a*ckg.wmk);  
  SNT snt_circ(m, r, circ_a, ckg.wmk, swsh_circ.lambda, Min(eps, 1e-6));
  CKR ckr(&swsh_circ, &snt_circ, &ckg, eps);
  Real EInfCirc =  abs(ckr.ZedH)*abs(ckr.ZedH) / (4.0*M_PI*iekr.omega*iekr.omega);

  // generic orbit 
  IEKG iekgINC(prograde,a,e,p,theta_,Min(eps, 1e-12));
  IEKR iekrINC(&iekgINC,l,m-k,k,n,eps);  // TEST  m --> m-k
  Complex ZHINC = iekrINC.ZH;
  Complex ZInfINC = iekrINC.ZInf;
  Real alphaINC = iekrINC.alpha;
  Real omegaINC = iekgINC.Omega(m-k,k,n); // TEST m --> m-k
  Real EInfINC = ZHINC.real()*ZHINC.real() + ZHINC.imag()*ZHINC.imag();
  EInfINC /= 4.0*M_PI*omegaINC*omegaINC;
  Real EHINC = ZInfINC.real()*ZInfINC.real() + ZInfINC.imag()*ZInfINC.imag();
  EHINC *= alphaINC;
  EHINC /= 4.0*M_PI*omegaINC*omegaINC;

  // inclined circular orbit
  CKG ckgINC(m-k, k, iekgINC.CosIota, USING_cosiota, r, circ_a); // TEST m --> m-k
  SWSH swsh_circINC(l, -2, m-k, circ_a*ckgINC.wmk);  // TEST m --> m-k
  SNT snt_circINC(m-k, r, circ_a, ckgINC.wmk, swsh_circINC.lambda, Min(eps, 1e-6)); // TEST m --> m-k
  CKR ckrINC(&swsh_circINC, &snt_circINC, &ckgINC, eps);
  Real EInfCircINC =  abs(ckrINC.ZedH)*abs(ckrINC.ZedH) / (4.0*M_PI*iekrINC.omega*iekrINC.omega);

  // compute Wigner-D function
  int mprime = m - k;   // use this if using wigner_d_www or wigner_d_arfken
  int mprime_old = k - m; // use this if using wigner_d_old
  //Real D = wigner_d_arfken(l, m, mprime, theta_inc);
  //Real D = wigner_d_old(l,m,mprime_old, theta_inc);
  //Real D = wigner_d_old(l,mprime, m, theta_inc);
  Real D = wigner_d_arfken(l,m,mprime, theta_inc);

  // compute residual
  circular_resid = fabs(   ( EInfCircINC / (EInfCirc*D*D) ) - 1.0   );
  generic_resid = fabs(   ( EInfINC / (EInf*D*D) ) - 1.0   );

  // TEST
  // display residuals
  /*
  cout.precision(10);
  Real residInf = fabs(   ( EInfINC / (EInf*D*D) ) - 1.0   );   
  Real residH = fabs(   ( EHINC / (EH*D*D) ) - 1.0   );
  cout << "%       EInfCircINC = " << EInfCircINC << "\t   EInfCirc = " << EInfCirc << endl;
  cout << "%           EInfINC = " << EInfINC     << "\t       EInf = " << EInf   << "\t D^2 = " << D*D << endl;
  //cout << "%         EHINC = " << EHINC       << "\t         EH = " << EH     << endl;
  cout << "%  generic_residual = " << residInf << endl; //  << "\t residual_H = " << residH << endl; 
  cout << "% circular_residual = " << fabs(   ( EInfCircINC / (EInfCirc*D*D) ) - 1.0   ) << endl;
  cout << "%     Scott's check = " << fabs(   ( EInfCircINC / EInfCirc ) - D*D   ) << endl;
  */

  // TEST
  // dump [theta_inc D*D]
  /*
  Real D_old, iota_rad;
  for(Real iota_deg = 1e-6; iota_deg <= 90.0; iota_deg += 90.0/500.0){
    iota_rad = iota_deg * M_PI / 180.0;
    D_old = wigner_d_old(l, m, mprime_old, iota_rad);
    D = wigner_d_arfken(l, m, mprime, iota_rad);
    cout << iota_deg << "\t" << D*D << "\t" << D_old*D_old << endl; 
  }
  */

  // normal exit
  return;

}

int main(int argc, char **argv)
{
 ios::sync_with_stdio();

  if(argc != 8) {
    cerr << "Arguments: 1. e  3. p " << endl;
    cerr << "           4. l  5. m  6. k  7. n" << endl;
    cerr << "           8. eps" << endl;
    exit(1);
  }
  const Real a = 0.0;
  const Real e = (Real)atof(argv[1]);
  const Real p = (Real)atof(argv[2]);
  int l = atoi(argv[3]);
  int m = atoi(argv[4]);
  int k = atoi(argv[5]);
  int n = atoi(argv[6]);
  Real eps = (Real)atof(argv[7]);

  // dump [theta_inc residual]
  cout.precision(10);
  Real generic_resid, circular_resid;
  for(Real theta_inc = 1e-6; theta_inc <= 90.0; theta_inc += 90.0/10.0){
    GetResidual(e, p, theta_inc, l, m, k, n, eps, circular_resid, generic_resid);
    cout << theta_inc << "\t" << circular_resid << "\t" << generic_resid << endl; 
  }
}
