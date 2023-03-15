//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Rec.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// This code will perform interpolations along a radiation reaction
// sequence to generate a trajectory through the parameter space of
// circular, equatorial Kerr geodesic orbits.
//
// Input parameters: just needs masses of the bodies and the filename.
// The initial radius is in the file, and the masses then set all the
// timescales.
//
// In particular, note that the time used in the integrator is related
// to Boyer-Lindquist time by t_BL = (M*M/mu) t_Code.  We actually
// output everything per unit mass, so we output \hat t_BL = (M/mu)
// t_Code = t_Code/eta.  In the final computation and output, this
// scaling is correctly accounted for in all quantities.
//
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include "Globals.h"
#include "CEID.h"
#include "NRUtil.h"

#define CLIGHT 2.9979248e5 // km/sec --- used to set the scale of the recoil.

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 9) {
    cerr << "Arguments: 1. data file basename" << endl;
    cerr << "           2. rmin  3. rmax  4. dr  5. rstart" << endl;
    cerr << "           6. lmax  7. M  8. mu" << endl;
    exit(1);
  }
  char basename[40], recname[40];
  FILE *infile, *recfile;
  //
  sprintf(basename, "%s", argv[1]);
  const Real rmin = (Real)atof(argv[2]);
  const Real rmax = (Real)atof(argv[3]);
  const Real dr = (Real)atof(argv[4]);
  const Real rstart = (Real)atof(argv[5]);
  const int lmax = atoi(argv[6]);
  const Real M = (Real)atof(argv[7]);
  const Real mu = (Real)atof(argv[8]);
  const Real eta = mu/M;
  //
  CEID ceid(basename, rmin, rmax, dr, lmax);
  ceid.Pfluxes();
  //
  // Start spiralling on in...
  //
  sprintf(recname, "%s.rec", basename);
  if (!(recfile = fopen(recname, "w"))) {
    cerr << "Error opening " << recname << endl;
    exit(2);
  }
  // bookkeepers
  int DONE_YET = 0, n = 0;
  //
  // dt is chosen so that we get good resolution on the inspiral.
  //
  Real dt = eta/(1000.*ceid.Omega_max);
  //
  // Initial conditions
  //
  Real t = 0., phi = 0., tau = 0.;
  Real rad = rstart;
  Real f_eta = eta*eta*sqrt(1. - 4.*eta);
  Real vx_cm, vy_cm, vz_cm;
  Real x_cm, y_cm, z_cm;
  fprintf(recfile, "#f(eta) = %e\n", f_eta);
  while (!DONE_YET) {
    ceid.Get_rdot_omega(rad);
    ceid.Get_fluxes(rad);
    ceid.Get_mom(rad, phi);
    if (n == 0) {
      vx_cm = ceid.pdot_y.real()*CLIGHT/ceid.Omega;
      vy_cm = -ceid.pdot_x.real()*CLIGHT/ceid.Omega;
      vz_cm = 0.;
      x_cm = -ceid.pdot_x.real()/(ceid.Omega*ceid.Omega);
      y_cm = -ceid.pdot_y.real()/(ceid.Omega*ceid.Omega);
      z_cm = 0.;
    }
    Real pdot_mag = sqrt(ceid.pdot_x.real()*ceid.pdot_x.real() +
			 ceid.pdot_y.real()*ceid.pdot_y.real());
    //
    //  1. (BL time)/M
    //  2. (radius)/M
    //  3. phi
    //  4. (radius)cos(phi)/M
    //  5. (radius)sin(phi)/M
    //  6. Edot
    //  7. Ldot
    //  8. rdot
    //  9. vx_cm/f
    // 10. vy_cm/f
    // 11. sqrt(vx_cm^2 + vy_cm^2)/f
    // 12. vx_cm
    // 13. vy_cm
    // 14. sqrt(vx_cm^2 + vy_cm^2)
    // 15. vfitchett
    // 16. x_cm/f
    // 17. y_cm/f
    // 18. sqrt(x_cm^2 + y_cm^2)/f
    // 19. x_cm
    // 20. y_cm
    // 21. sqrt(x_cm^2 + y_cm^2)
    // 22. sqrt(pdot_x^2 + pdot_y^2)/f
    // 23. sqrt(pdot_x^2 + pdot_y^2)
    //
    if (n%50 == 0)
      fprintf(recfile,
	      "%e %lf %lf %lf %lf %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	      t, rad, phi, rad*cos(phi), rad*sin(phi),
	      ceid.Edot*eta, ceid.Lzdot*eta, ceid.rdot*eta,
	      vx_cm, vy_cm, sqrt(vx_cm*vx_cm + vy_cm*vy_cm),
	      f_eta*vx_cm, f_eta*vy_cm, f_eta*sqrt(vx_cm*vx_cm + vy_cm*vy_cm),
	      CLIGHT*29.*f_eta*pow(2./rad, 4.)/105.,
	      x_cm, y_cm, sqrt(x_cm*x_cm + y_cm*y_cm),
	      f_eta*x_cm, f_eta*y_cm, f_eta*sqrt(x_cm*x_cm + y_cm*y_cm),
	      pdot_mag, f_eta*pdot_mag);
    //
    // Ending condition: If dt gets way tiny, kill the loop.  We
    // decrease dt if the next step will go inside the minimum
    // radius.  This radius is r[imax].
    //
    if (dt > 1.e-7) {
      Real dr_tmp = ceid.rdot*dt;
      while (rad + dr_tmp < ceid.rmin && dt > 1.e-8) {
	dt /= 2.;
	dr_tmp = ceid.rdot*dt;
      }
      //
      // Once through there, update the coordinates.
      //
      t += dt/eta;
      phi += ceid.Omega*dt/eta;
      rad += dr_tmp;
      vx_cm += ceid.pdot_x.real()*CLIGHT*dt/eta;
      vy_cm += ceid.pdot_y.real()*CLIGHT*dt/eta;
      vz_cm += ceid.pdot_z.real()*CLIGHT*dt/eta;
      x_cm += vx_cm*dt/(CLIGHT*eta);
      y_cm += vy_cm*dt/(CLIGHT*eta);
      z_cm += vz_cm*dt/(CLIGHT*eta);
      n++;
    } else {
      DONE_YET = 1;
    }
  }
  fclose(recfile);
}

