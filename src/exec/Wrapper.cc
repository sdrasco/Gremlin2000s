//---------------------------------------------------------------------------
//
// $Id: Wrapper.cc,v 1.27 2004/09/04 19:28:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include "Globals.h"
#include "IEKG.h"
#include "NRIEKG.h"
#include "CKG.h"
#include "SNT.h"
#include "SWSH.h"
#include "CKR.h"
#include "IEKR.h"

void compare(Real p,Real e,Real thetaminus,Real a,
             int l,int m,int k,int n,Real eps,
             Real &CosIota, Real &FracErr)
{

  //
  // Create Inclined Eccentric Kerr Geodesic
  //
  IEKG iekg(1, a, e, p, thetaminus,eps);
  Real iekgwmkn = iekg.Omega(m,k,n);
  //
  // Chose an good radius for circular orbits to compare with
  //
  Real MinR = p/(1.0 + e); 
  Real MaxR = p/(1.0 - e);
  Real AverageR = 0.5*(MinR+MaxR);
  Real r = AverageR;
  //
  // Create Circular Kerr Geodesic
  //
  CKG ckg(m, k, iekg.CosIota, USING_cosiota, r, a);
  //
  // Create spheroidal harmonics for both orbital frequencies
  //
  SWSH swsh_ecc(l, -2, m, a*iekgwmkn);
  SWSH swsh_circ(l, -2, m, a*ckg.wmk);
  //
  // Homogeneous Sasaki-Nakamura-Teukolsky for both
  //
  SNT snt_ecc(m, r, a, iekgwmkn, swsh_ecc.lambda, Min(eps, 1e-6));
  SNT snt_circ(m, r, a, ckg.wmk, swsh_circ.lambda, Min(eps, 1e-6));
  //
  // Circular Kerr Radiation
  //
  CKR ckr(&swsh_circ, &snt_circ, &ckg, eps);
  //
  // Inclined Eccentric Kerr Radiation
  //
  IEKR iekr(&iekg,l,m,k,n,eps);
  //
  // Return our standard comparison, and inclination angle
  //
  FracErr = abs(iekr.ZH - ckr.ZedH) / abs(iekr.ZH);
  CosIota = iekg.CosIota;
  //
  // optionally spit out some other comparison data
  //
  Real EdotInfGeneric = abs(iekr.ZH)*abs(iekr.ZH) / (4.0*M_PI*iekr.omega*iekr.omega);
  Real EdotInfCircular =  abs(ckr.ZedH)*abs(ckr.ZedH) / (4.0*M_PI*iekr.omega*iekr.omega);
  FracErr = fabs(EdotInfGeneric-EdotInfCircular) / EdotInfGeneric;
  cout.precision(16);
  /*
  cout << k << "\t"
       << fabs(EdotInfGeneric-EdotInfCircular) << "\t"
       << EdotInfGeneric << "\t"
       << endl;
  */
  cout << "\n\tFractional error is: " << FracErr << "\n\n";
  /*
  cout << endl
       //<< "Steve: omega = " << iekgwmkn << endl
       //<< "Scott: omega = " << ckg.wmk << endl
       //<< "Steve: TeukRH = " << snt_ecc.TeukRH << endl
       //<< "Scott: TeukRH = " << snt_circ.TeukRH << endl
       //<< "Steve: TeukRInf = " << snt_ecc.TeukRInf << endl
       //<< "Scott: TeukRInf = " << snt_circ.TeukRInf << endl
       << "Steve: ZH = " << iekr.ZH << endl
       << "Scott: ZH = " << ckr.ZedH << endl
       << "Steve: ZInf = " << iekr.ZInf << endl
       << "Scott: ZInf = " << ckr.ZedInf << endl
       //<< "theta_ = " << thetaminus << "\t" 
       //<< "DeltaZH.real / ZH.real" << fabs(iekr.ZH.real() - ckr.ZedH.real()) / fabs(iekr.ZH.real()) 
       << endl;
  */

}

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  char outname[50];
  //
  if(argc != 10) {
    cerr << "Arguments: 1. p  2. e  3. Steve's inclination parameter (in degrees)" << endl;
    cerr << "           4. a  5. l  6. m  6. k  8. n" << endl;
    cerr << "           9. eps" << endl;
    cerr << "inclined-circular example: "
         << "bin/Wrapper 7 1e-7 45.0 0.998 2 2 0 0 1e-5" << endl;
    cerr << "equatorial-circular example: "
         << "bin/Wrapper 7 1e-7 89.99999 0.998 2 2 0 0 1e-5" << endl;
    exit(1);
  }
  const Real p = (Real)atof(argv[1]);
  const Real e = (Real)atof(argv[2]);
  const Real thetaminus_deg = (Real)atof(argv[3]);
  const Real thetaminus = thetaminus_deg * M_PI/180.0;
  const Real a = (Real)atof(argv[4]);
  const int l = atoi(argv[5]);
  const int m = atoi(argv[6]);
  const int k = atoi(argv[7]);
  const int n = atoi(argv[8]);
  const Real eps = (Real)atof(argv[9]);

  // get IE and CIRC comparison figure for the requested inclination
  Real CosIota;
  Real FracErr;
  compare(p,e,thetaminus,a,l,m,k,n,eps,CosIota,FracErr);
  /*
  cout.precision(16);
  cout << acos(CosIota) << "\t"
       << FracErr << "\t"
       << "\n";
  */

  //
  // forget the requested inclination, do a bunch of them (put in by hand)
  //
  // lots of inclinations
  /*
  Real CosIota;
  Real FracErr;
  Real start = 0.1;
  Real finish = M_PI/2.0 - eps;
  Real samples = 500;
  Real deltax = (finish - start) / samples;
  for(Real x = start; x < finish; x += deltax) {
    compare(p,e,x,a,l,m,k,n,eps,CosIota,FracErr);
    cout.precision(16);
    cout << acos(CosIota) << "\t"
         << FracErr << "\t"
         << x << "\t"
         << "\n";
    cout.flush();
  }
  */
  //
  // lots of lmkn
  /*
  Real CosIota;
  Real FracErr;
  cout << "% k  |       abs. resid.      |    Edot\n";
  for(int kk=-10; kk <= 10; kk++){
    compare(p,e,thetaminus,a,l,m,kk,n,eps,CosIota,FracErr);
  }
  */
}
  

