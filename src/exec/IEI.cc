//   Test code for IETD: evolution of generic orbits
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

void GetXdot(Real a, Real &e, Real &p, Real &theta_inc, Real E, Real Lz, Real Q, Real &dEH, Real &dEI, Real &dLH, Real &dLI, Real &dQ, Real &dt, Real eps)
{

  IEKG iekg(a, E, Lz, Q, 1e-6);
  e = iekg.e;
  p = iekg.p;
  if(iekg.prograde) {
    theta_inc = 90.0 - (180.0/M_PI)*iekg.theta_;
  } else {
    theta_inc = 90.0 + (180.0/M_PI)*iekg.theta_;
  }
  dt = (500.0)*2.0*M_PI / fabs(iekg.OmegaPhi);

  IETD ietd(&iekg, eps);
  dEH = ietd.dEH;
  dEI = ietd.dEI;
  dLH = ietd.dLH;
  dLI = ietd.dLI;
  dQ = ietd.dQ;
}

int main(int argc, char **argv)
{

  ios::sync_with_stdio();
  char outname[50];

  if(argc != 8) {
    cerr << "Arguments: 1. M/Msun  2. mu/Msun 3. a  4. e_initial  5. p_initial " << endl;
    cerr << "           6. theta_min,initial (in degrees) " << endl;
    cerr << "           7. eps" << endl;
    exit(1);
  }

  // read inputs
  Real M = (Real)atof(argv[1]);
  Real mu = (Real)atof(argv[2]);
  Real a = (Real)atof(argv[3]);
  Real e = (Real)atof(argv[4]);
  Real p = (Real)atof(argv[5]);
  Real MinTheta = (Real)atof(argv[6]);
  Real eps = (Real)atof(argv[7]);
  Real theta_ = M_PI * fabs(MinTheta)/180.0;
  int prograde = (MinTheta > 0 ? 1 : 0);

  // compute initial constants of motion
  Real eta = mu/M;
  IEKG iekg(prograde,a,e,p,theta_,Min(eps, 1e-12));
  Real E = iekg.E;
  Real L = iekg.Lz;
  Real Q = iekg.Q;
  Real dEH, dEI, dLH, dLI, dQ;

  // compute Edot_isco (see Finn & Thorne).
  // NOTE: this is the normalized version of Edot_isco shown in Finn & Thorne.
  Real Edot_isco = (32.0/5.0) * pow(iekg.Omega_isco,10.0/3.0);
  
  // compute timestep from Edot_isco
  //Real dt = 0.5*(M/mu) * E / Edot_isco;

  // write header
  cout << "%\n";
  cout << "% a/M = " << iekg.a << "\n";
  cout << "% M/Msun = " << M << "\n";
  cout << "% mu/Msun = " << mu << "\n";
  cout << "% r_isco/M = " << iekg.r_isco << "\n";
  cout << "% Edot_isco = " << Edot_isco << "\n";
  //cout << "% dt/M = " << dt << "\n";
  cout << "%\n";
  cout << "% output format [t e p theta_inc (deg) E Lz Q <dE/dt>_H <dE/dt>_Inf <dLz/dt>_H <dLz/dt>_Inf <dQ/dt>] \n";
  cout << "%\n";
  cout << "% NOTE: This table shows <dX/dt> = - <dX/dt>_particle.  Even in the case of <dQ/dt>.\n";
  cout << "%\n";

  // TEST try using orbital timescale to set timestep
  Real t = 0.0;
  Real dt;
  
  // main loop
  Real theta_inc;
  //for(int i=1; i<=10; i++){
  int i = 0;
  Real p_SchSep = 6.0 + 2*e;
  //while(p >  1.1*p_SchSep){
  while(p>2){

    // increment index
    i++;

    // radiate
    GetXdot(a, e, p, theta_inc, E, L, Q, dEH, dEI, dLH, dLI, dQ, dt, eps);
    cout.precision(16);
    //cout << (Real(i)-1.0)*dt << "\t" 
    cout << t+dt << "\t"
         << e << "\t"
         << p << "\t"
         << theta_inc << "\t"
         << E << "\t"
         << L << "\t"
         << Q << "\t"
         << dEH << "\t"
         << dEI << "\t"
         << dLH << "\t"
         << dLI << "\t"
         << dQ << "\t"
         << "\n";
    cout.flush();

    // compute timestep from initial energy flux
    /*
    if(i == 1) {
      Real Edot_initial = dEH + dEI;  
      Real dt = (1.0e-7) * (M/mu) * E / Edot_initial;
    }
    */
    t += dt;

    // evolve
    E = E - dt*(mu/M)*(dEH + dEI);
    L = L - dt*(mu/M)*(dLH + dLI);
    Q = Q - dt*(mu/M)*dQ;

  }

}
