//   For simple tests of the IEKG stuff
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
#include "IETD.h"

// constant
Real TEPS = 0.001;

// little piece of code that returns R^H(r)
Real GetRH(IEKG *iekg, SWSH *swsh, Real r){
  int m = 2;
  int k = 1;
  int n = 1;
  SNT snt(m, r, iekg->a, iekg->Omega(m,k,n), swsh->lambda, 1e-6);
  return snt.TeukRH.real();
}


// little piece of code to return frequencies
void Freqs(int prograde,Real a,Real e,Real p,Real theta_,Real eps,
          Real *gamma, Real *OmegaR, Real *OmegaTheta, Real *OmegaPhi){

  // make the geodesic object
  IEKG iekg(prograde,a,e,p,theta_,eps);

  // assign output
  *gamma = iekg.gamma;
  *OmegaR = iekg.OmegaR;
  *OmegaTheta = iekg.OmegaTheta;
  *OmegaPhi = iekg.OmegaPhi;

  // done
  return;
}

// compute waveform and flux information (H only for now)
void GetXlmkn(IEKG *iekg, int l, int m, int k, int n, Real eps)
{
  // (m,k,n) = (0,0,0) is a degenerate case
  //if(m==0 && k==0 && n==0) {
  if(fabs(iekg->Omega(m,k,n)) <= 1e-10){
    cout << "omega = " << iekg->Omega(m,k,n) << " is too small. All fluxes vanish\n";
    return;
  } else {
    // compute radiation
    IEKR iekr(iekg,l,m,k,n,eps);
    Complex ZH = iekr.ZH;
    Complex ZInf = iekr.ZInf;
    Real alpha = iekr.alpha;
    Real omega = iekg->Omega(m,k,n);

    // compute spheroidal harmonic (at a few different angles)
    SWSH swsh(l, -2, m, (iekg->a)*omega);
    Real ct90 = cos(M_PI/2.0 - TEPS);
    Real ct60 = cos(M_PI/3.0);
    Real ct30 = cos(M_PI/6.0);
    Real ct0 = cos(TEPS);
    Real S90 = swsh.spheroid(ct90);
    Real S60 = swsh.spheroid(ct60);
    Real S30 = swsh.spheroid(ct30);
    Real S0 = swsh.spheroid(ct0);

    // compute EInf
    Real EInf = ZH.real()*ZH.real() + ZH.imag()*ZH.imag();
    EInf /= 4.0*M_PI*omega*omega;

    // compute EH
    Real EH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
    EH *= alpha;
    EH /= 4.0*M_PI*omega*omega;

    // compute LInf
    Real LInf = ZH.real()*ZH.real() + ZH.imag()*ZH.imag();
    LInf *= Real(m);
    LInf /= 4.0*M_PI*omega*omega*omega;

    // compute LH
    Real LH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
    LH *= Real(m);
    LH *= alpha;
    LH /= 4.0*M_PI*omega*omega*omega;

    // record radiation data
    cout.precision(6);
    cout << "% "
         << l << "\t"
         << m << "\t"
         << k << "\t"
         << n << "\t"
         << omega << "\t"
         //<< ZH.real() << "\t"
         //<< ZH.imag() << "\t"
         //<< ZInf.real() <<  "\t"
         //<< ZInf.imag() << "\t"
         //<< iekr.alpha << "\t"
         //<< S90 << "\t"
         //<< S60 << "\t"
         //<< S30 << "\t"
         //<< S0 << "\t"
         //<< "|ZH|^2 = " << abs(ZH)*abs(ZH) << "\t"
         //<< "|ZI|^2 = " << abs(ZInf)*abs(ZInf) << "\t"
         //<< "alpha = " << alpha << "\t"
         << "EI = " << EInf << "\t"
         << "EH = " << EH << "\t"
         << "LI = " << LInf << "\t"
         << "LH = " << LH << "\t"
         << endl;
    cout.flush();
    }
}


int main(int argc, char **argv)
{

  ios::sync_with_stdio();
  char outname[50];
  /*
  //
  if(argc != 2) {
    cerr << "Arguments: 1. flag: 1 = equatorial\n"  
         << "                    2 = circular\n"
         << "                    3 = inclined and eccentric\n";
    exit(1);
  }
  const int flag = atoi(argv[1]);
  */
  //

  if(argc != 10) {
    cerr << "Arguments: 1. a  2. e  3. p " << endl;
    cerr << "           4. theta_min (in degrees)  5. l  6. m  6. k  8. n" << endl;
    cerr << "           9. eps" << endl;
    exit(1);
  }
  const Real a = (Real)atof(argv[1]);
  const Real e = (Real)atof(argv[2]);
  const Real p = (Real)atof(argv[3]);
  const Real MinTheta = (Real)atof(argv[4]);
  const int l = atoi(argv[5]);
  const int m = atoi(argv[6]);
  const int k = atoi(argv[7]);
  const int n = atoi(argv[8]);
  const Real eps = (Real)atof(argv[9]);

  // try out IEKG
  //const Real MinTheta = atof(argv[1]);
  const int prograde=1;
  //const Real a = 0.998, e = 3e-5, p = 2.0,eps=1e-5;
  const Real theta_ = MinTheta *M_PI / 180.0;
  IEKG iekg(prograde,a,e,p,theta_,eps);
  //cout << "\n\t Delta Q / Q = " << (iekg.Q - 1.316861192627782) / 1.316861192627782 << endl;

  // dump gamma Omega_r Omega_theta Omega_phi for a plot
  // edit IEKG before running this bit.  we don't want the
  // full orbits, just the frequencies.
  /*
  int prograde = 0;
  Real a = 0.9, e = 0.7, theta_ = 10*M_PI/180.0, eps=1e-6;
  //Real a = 0.998, e = 1.0/3.0, theta_ = M_PI/3.0, eps=1e-6;
  cout << "%a = " << a << "\n"
       << "%e = " << e << "\n"
       << "%theta_ = " << theta_ << "\n"
       << "%eps =" << eps << "\n"
       << "%[Rmin/M gamma M*OmegaR/(2pi) M*OmegaTheta/(2pi) M*OmegaPhi/(2pi) p]\n";
  Real start = 13;
  Real finish = 5.5;
  Real gamma, OmegaR, OmegaTheta, OmegaPhi, Rmin;
  //Real Delta = (log(start) - log(finish)) / 1000.0; // log spacingi
  //for(Real p = start; p >= finish; p *= exp(-Delta) ){
  Real Delta = (start - finish) / 100.0;
  for(Real p = start; p >= finish; p -= Delta ){
    Rmin = p / (1.0 + e);
    Freqs(prograde,a,e,p,theta_,eps,&gamma,&OmegaR,&OmegaTheta,&OmegaPhi);
    cout.precision(16);
    cout << Rmin << "\t"
         << gamma   << "\t"
         << OmegaR/(2.0*M_PI)  << "\t"
         << OmegaTheta/(2.0*M_PI) << "\t"
         << OmegaPhi/(2.0*M_PI) << "\t"
         << p << "\t"
         << endl;
    cout.flush();
  }
  */

  // dump simple figures 
  //
  // 20 June 2003: passes comparison to Wolfram Schmidt's published figures
  // 11 Nov 2003: passes comparison to Scott's astro-ph/0308479 code.
  // 
  cout.precision(16);
  cout << "%-------------------Parameter List-------------------\n";
  cout << "%IEKG.prograde = " << iekg.prograde << endl;
  cout << "%IEKG.a = " << iekg.a << endl;
  cout << "%IEKG.e = " << iekg.e << endl; 
  cout << "%IEKG.p = " << iekg.p << endl;
  cout << "%IEKG.theta_ = " << iekg.theta_ << endl;
  cout << "%IEKG.theta_ = " << MinTheta << " (degrees) " << endl;
  cout << "%IEKG.CosIota = " << iekg.CosIota << endl;
  cout << "%Iota = " << acos(iekg.CosIota) << endl;
  cout << "%IEKG.E = " << iekg.E << endl;
  cout << "%IEKG.Lz = " << iekg.Lz << endl;
  cout << "%IEKG.Q = " << iekg.Q << endl;
  cout << "%IEKG.OmegaR = " << iekg.OmegaR << endl;
  cout << "%IEKG.OmegaTheta = " << iekg.OmegaTheta << endl;
  cout << "%IEKG.OmegaPhi = " << iekg.OmegaPhi << endl;
  cout << "%IEKG.gamma = " << iekg.Gamma << endl;
  cout << "%IEKG.UpsilonR = " << iekg.UpsilonR << endl;
  cout << "%IEKG.UpsilonTheta = " << iekg.UpsilonTheta << endl;
  cout << "%IEKG.UpsilonPhi = " << iekg.UpsilonPhi << endl;
  cout << "%---------------------------------------------------\n";

  // dump w_r(psi) and w_theta(chi) as a test
  // 
  // 11 Nov 2003: passes comparison to Scott's astro-ph/0308479 code.
  // 
  /*
  Real x;
  int N = 500;
  cout.precision(16);
  for(int i=0; i<N; i++){
    x = 2.0*M_PI* double(i)/double(N);
    cout << x << "\t";
    cout << iekg.wR(x) << "\t";
    cout << iekg.wTheta(x) << endl;
  }
  */

  // dump phi(lambda) and t(lambda)   
  //
  // 17 Feb 2004: passes comparison to Scott's public geodesic code.
  // 29 Feb 2004: shift t(lambda) by -0.08 when comparing with Scott's
  //              "Generic" routine.
  // 
  cout.precision(16);
  Real t = 0.0, phi = 0.0;
  Real LambdaMax = 3.0 * 2.0*M_PI / iekg.UpsilonR;
  Real DeltaLambda = LambdaMax / 10000.0;
  cout << "% dumping below: [lambda*M t/M phi(radians)]" << endl << "% " << endl;
  for(Real Lambda = 0.0; Lambda < LambdaMax; Lambda += DeltaLambda){
    t = 0.0;
    phi = 0.0;
    for(int k=1; k <= iekg.Delta_t_kmax; k++) {
      t+=iekg.Delta_t_k[k]*sin(double(k)*iekg.UpsilonTheta*Lambda);
    }
    for(int n=1; n <= iekg.Delta_t_nmax; n++) {
      t+=iekg.Delta_t_n[n]*sin(double(n)*iekg.UpsilonR*Lambda);
    }
    for(int k=1; k <= iekg.Delta_phi_kmax; k++) {
      phi+=iekg.Delta_phi_k[k]*sin(double(k)*iekg.UpsilonTheta*Lambda);
    }
    for(int n=1; n <= iekg.Delta_phi_nmax; n++) {
      phi+=iekg.Delta_phi_n[n]*sin(double(n)*iekg.UpsilonR*Lambda);
    }
    t += iekg.Gamma * Lambda;
    phi += iekg.UpsilonPhi * Lambda;
    cout << Lambda << "\t" << t << "\t" << phi << endl;
  } 

  // dump r(lambda)
  /*
  cout.precision(16);
  Real PsiMax = 0.9*M_PI;
  Real DeltaPsi = PsiMax / 1000.0;
  for(Real Psi=0.0; Psi <= PsiMax; Psi += DeltaPsi)
    cout << wR(Psi,iekg.args)/iekg.UpsilonR << "\t" << RofPsi(Psi,iekg.args) << "\n";
  */

  // dump theta(lambda)
  /*
  cout.precision(16);
  Real ChiMax = 0.9*M_PI;
  Real DeltaChi = ChiMax / 1000.0;
  for(Real Chi=0.0; Chi <= ChiMax; Chi += DeltaChi) 
    cout << iekg.wTheta(Chi)/iekg.UpsilonTheta 
         << "\t" << -cos(iekg.Theta(Chi)) << "\n";
    // the -cos(theta) makes my output have the same initial condition as Generic.
  */

  // dump Delta_t to check symmetry
  /*
  cout.precision(16);
  Real XMax = 2.0*M_PI;
  Real DeltaX = XMax / 10000.0;
  for(Real X = 0.0; X < XMax; X += DeltaX){
    cout << X << "\t" << iekg.Delta_t(1.0,X) << endl;
  }
  */


  // try out IEKR 
  /*
  Real theta_ = M_PI * fabs(MinTheta)/180.0;
  int prograde = (MinTheta > 0 ? 1 : 0);
  IEKG iekg(prograde,a,e,p,theta_,Min(eps, 1e-12));
  GetXlmkn(&iekg, l, m, k, n, eps);
  */

  // try overloaded IEKG(a,E,Lz,Q,eps)  (only input argument is a)
  /*
  Real E =  0.906207487335045;
  Real L = 2.532133117162242;
  Real Q = 2.470768634192265e-32;
  IEKG iekg(a, E, L, Q, 1e-12);
  cout.precision(16);
  cout << "e = " << iekg.e << "\n" 
       << "p = " << iekg.p << "\n"
       << "theta_ = " << (180.0/M_PI)*iekg.theta_ << "\n";
  */

  // test of overloaded IEKG(a,E,Lz,Q,eps)
  /*
  Real theta_ = M_PI * fabs(MinTheta)/180.0;
  int prograde = (MinTheta > 0 ? 1 : 0);
  IEKG iekg(prograde,a,e,p,theta_,eps);
  cout << "<cot^2 theta> = " << iekg.avgcot2 << endl;
  */
  /*
  IEKG iekg2(a, iekg.E, iekg.Lz, iekg.Q, 1e-12);
  cout.precision(16);
  cout << "E = " << iekg.E << "\n" 
       << "Lz = " << iekg.Lz << "\n"
       << "Q = " << iekg.Q << "\n";
  cout << "e = " << iekg2.e << "\n" 
       << "p = " << iekg2.p << "\n"
       << "theta_ = " << (180.0/M_PI)*iekg2.theta_ << "\n";
  cout << "Delta e / e = " << (iekg.e - iekg2.e) / iekg.e << "\n"
      << "Delta theta_ / theta_ = " << (iekg.theta_ - iekg2.theta_) / iekg.e << "\n"
      << "Delta p / p = " << fabs(iekg.p - iekg2.p) / iekg.p << "\n"
      << "\n";
  cout.flush();
  GetXlmkn(&iekg, l, m, k, n, eps);
  GetXlmkn(&iekg2, l, m, k, n, eps);
  */

  // try out IETD
  /*
  Real theta_ = M_PI * fabs(MinTheta)/180.0;
  int prograde = (MinTheta > 0 ? 1 : 0);
  IEKG iekg(prograde,a,e,p,theta_,Min(eps, 1e-10));
  IETD ietd(&iekg, eps);
  cout.precision(16);
  cout << ietd.dEH << "\t"
      << ietd.dEI << "\t"
      << ietd.dLH << "\t"
      << ietd.dLI << "\t"
      << ietd.dQ << "\t"
      << "\n";
  cout.flush();
  */

  // practice with Complex
  //Complex x = Complex(0.0,M_PI);
  //Complex z = Complex(0.0,1.0);
  //cout << x.imag() << endl;
  //cout << x.real() << endl;
  //cout << exp(x) << endl;
  //Complex *y = &x;
  //*y = z;
  //cout << y[0] << endl; 
  //cout << x << endl;

  //
  //  Each of the following have nice doable orbits
  //  
  /*
  int klimit = 3;
  int nlimit = 3;
  if (flag == 3){
    //
    //  inclined and eccentric (do a loop over both k and n)
    //
    const int prograde=1;
    const Real a = 0.9, e = 0.5, MinTheta = 45.0001, p = 4.0, eps=1e-5;
    Real theta_ = MinTheta * M_PI/180.0;
    IEKG iekg(prograde,a,e,p,theta_,eps);
    cout.precision(16);
    cout << "%\n"
         << "%a = " << a << "\n"
         << "%e = " << e << "\n"
         << "%p = " << p << "\n"
         << "%MinTheta = " << MinTheta << " (degrees)\n"
         << "%eps = " << eps << "\n"
         << "%Teps = " << TEPS << "\n"
         << "%Delta_t_kmax = " << iekg.Delta_t_kmax << "\n"
         << "%Delta_t_nmax = " << iekg.Delta_t_nmax << "\n"
         << "%Delta_phi_kmax = " << iekg.Delta_phi_kmax << "\n"
         << "%Delta_phi_nmax = " << iekg.Delta_phi_nmax << "\n"
         << "%\n"
         << "%below: [l m k n omega Z^H.real() ZH.imag() Z^Inf.real() ZInf.imag() ... \n"
         << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH]" << "\n"
         << "%\n";
    cout.flush();
    int l = 2;
    int m = 2;
    Complex ZH;
    for(int k = -klimit; k <= klimit; k++){
      for(int n = -nlimit; n<= nlimit; n++){
        GetXlmkn(&iekg, l, m, k, n, eps);
      }
    }
  } else if (flag == 2) {
    //
    //  inclined, but not eccentric (do a loop over k)
    //  
    const int prograde=1;
    const Real a = 0.9, e = 1e-7, MinTheta = 45.0001, p = 4.0, eps=1e-5;
    Real theta_ = MinTheta * M_PI/180.0;
    IEKG iekg(prograde,a,e,p,theta_,eps);
    cout.precision(16);
    cout << "%\n"
         << "%a = " << a << "\n"
         << "%e = " << e << "\n"
         << "%p = " << p << "\n"
         << "%MinTheta = " << MinTheta << " (degrees)\n"
         << "%eps = " << eps << "\n"
         << "%Teps = " << TEPS << "\n"
         << "%Delta_t_kmax = " << iekg.Delta_t_kmax << "\n"
         << "%Delta_t_nmax = " << iekg.Delta_t_nmax << "\n"
         << "%Delta_phi_kmax = " << iekg.Delta_phi_kmax << "\n"
         << "%Delta_phi_nmax = " << iekg.Delta_phi_nmax << "\n"
         << "%\n"
         << "%below: [l m k n omega Z^H.real() ZH.imag() Z^Inf.real() ZInf.imag() ... \n"
         << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH]" << "\n"
         << "%\n";
    cout.flush();
    int l = 2;
    int m = 2;
    int n = 0;
    Complex ZH;
    for(int k = -klimit; k <= klimit; k++){
      GetXlmkn(&iekg, l, m, k, n, eps);
    }
  } else if (flag == 1) {
    //
    //  not inclined, but eccentric (do a loop over n)
    //
    //const int prograde=1;
    //const Real a = 0.9, e = 0.5, p = 4.0, eps=1e-5;
    // 
    // orbit from Glampedakis & Kennifick Figs. 14 and 15.
    //
    const int prograde=1;
    Real a = 0.99, e = 0.7, p = 2.11, eps=1e-8;
    //
    // orbit from Glampedakis & Kennifick Fig. 16.
    //
    //const int prograde=0;
    //Real a = 0.99, e = 0.5, p = 10.4, eps=1e-6;
    //
    //
    Real MinTheta = 90.0 - eps;
    Real theta_ = MinTheta * M_PI/180.0;
    IEKG iekg(prograde,a,e,p,theta_,eps);
    cout.precision(16);
    cout << "%\n"
         << "%a = " << a << "\n"
         << "%e = " << e << "\n"
         << "%p = " << p << "\n"
         << "%MinTheta = " << MinTheta << " (degrees)\n"
         << "%eps = " << eps << "\n"
         << "%Teps = " << TEPS << "\n"
         << "%Delta_t_kmax = " << iekg.Delta_t_kmax << "\n"
         << "%Delta_t_nmax = " << iekg.Delta_t_nmax << "\n"
         << "%Delta_phi_kmax = " << iekg.Delta_phi_kmax << "\n"
         << "%Delta_phi_nmax = " << iekg.Delta_phi_nmax << "\n"
         << "%\n"
         //<< "%below: [l m k n omega Z^H.real() ZH.imag() Z^Inf.real() ZInf.imag() ... \n"
         //<< "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH]" << "\n"
         << "%below: [l m n omega  EInf EH LInf LH]\n"
         << "%\n";
    cout.flush();
    int l = 2;
    int m = 2;
    int k = 0;
    int min_n = -20, max_n = 20;
    Complex ZH;
    for(int n = min_n; n<= max_n; n++){
      GetXlmkn(&iekg, l, m, k, n, eps);
    }
  } else { 
    cout << "Unknown flag: 1 = equatorial;  2 = circular; 3 = inclined and eccentric;\n";
  }
  */

  // test out new version of SNT
  /*
  const Real a = 0.99, e = 0.5, theta_ = 0.7853981, p = 4.0, eps=1e-5;
  int l=5,m=2,n=10,k=10,prograde=1;
  IEKG iekg(prograde,a,e,p,theta_,eps);
  Real iekgwmkn = iekg.Omega(m,k,n);
  SWSH swsh(l, -2, m, a*iekgwmkn);
  int NumSteps = 2000;
  SNT snt(m, iekg.p, iekg.e, a, iekgwmkn, swsh.lambda, Min(eps, 1e-6));
  Real rmin = p/(1.0+e);
  Real rmax = p/(1.0-e);
  Real r, dr = (rmax - rmin) / Real(NumSteps-1);
  cout.precision(16);
  for(int i=0; i<NumSteps; i++){
    r = rmin + Real(i)*dr;
    snt.CalcRXFields(r,1);
    cout << r << "\t" << snt.TeukRH.real() << endl;
  }
  int NumSpline = 20000;
  snt.PrepSpline(NumSpline);
  for(int i=0; i<NumSteps; i++){
    r = rmin + Real(i)*dr;
    snt.SplineFields(r,1);
    cout << r << "\t" << snt.TeukRH.real() << endl;
  }
  */

  // get R^H the hard way (for comparison to new SNT)
  /*
  const int prograde=1;
  const Real a = 0.99, e = 0.5, theta_ = 0.78540, p = 4.0, eps=1e-5;
  int l=2,m=2,n=1,k=1;
  IEKG iekg(prograde,a,e,p,theta_,eps);
  Real iekgwmkn = iekg.Omega(m,k,n);
  SWSH swsh(l, -2, m, a*iekgwmkn);
  int NumSteps = 2000;
  Real rp = 1.0 + sqrt(1.0 - a*a);
  Real dr = (iekg.rmax - rp) / Real(NumSteps);
  cout.precision(16);
  for(Real r=rp; r<iekg.rmax; r+=dr){
    cout << r << "\t" << GetRH(&iekg, &swsh, r) << endl;
    cout.flush();
  }
  */

}
