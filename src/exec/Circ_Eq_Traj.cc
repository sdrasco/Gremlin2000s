//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Traj.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
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

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 8) {
    cerr << "Arguments:  1. data file basename" << endl;
    cerr << "            2. rmin  3. rmax  4. dr" << endl;
    cerr << "            5. rstart  6. M  7. mu" << endl;
    exit(1);
  }
  char basename[40], trajname[40];
  FILE *trajfile;
  //
  sprintf(basename, "%s", argv[1]);
  const Real rmin = (Real)atof(argv[2]);
  const Real rmax = (Real)atof(argv[3]);
  const Real dr = (Real)atof(argv[4]);
  const Real rstart = (Real)atof(argv[5]);
  const Real M = (Real)atof(argv[6]);
  const Real mu = (Real)atof(argv[7]);
  const Real eta = mu/M;
  //
  CEID ceid(basename, rmin, rmax, dr, 0);
  //
  // Start spiralling on in...
  //
  sprintf(trajname, "%s_r%2.1lf.traj", basename, rstart);
  if (!(trajfile = fopen(trajname, "w"))) {
    cerr << "Error opening " << trajname << endl;
    exit(3);
  }
  // bookkeeper
  int DONE_YET = 0;
  //
  // dt is chosen so that we get good resolution on the inspiral.
  //
  Real dt = eta/(10.*ceid.Omega_max);
  //
  // Initial conditions
  //
  Real t = 0., phi = 0.;
  Real rad = rstart;
  while (!DONE_YET) {
    ceid.Get_rdot_omega(rad);
    //
    //  1. (BL time)/M
    //  2. (radius)/M
    //  3. phi
    //  4. (radius)cos(phi)/M
    //  5. (radius)sin(phi)/M
    //
    fprintf(trajfile, "%e %lf %lf %lf %lf\n",
	    t, rad, phi, rad*cos(phi), rad*sin(phi));
    //
    // Ending condition: If dt gets way tiny, kill the loop.  We
    // decrease dt if the next step will go inside the minimum
    // radius.
    //
    if (dt > 1.e-7) {
      Real dr_tmp = ceid.rdot*dt;
      while (rad + dr_tmp < ceid.rmin && dt > 1.e-8) {
	dt /= 2.;
	dr_tmp = ceid.rdot*dt;
      }
      //
      // Update the coordinates.
      //
      t += dt/eta;
      phi += ceid.Omega*dt/eta;
      rad += dr_tmp;
    } else {
      DONE_YET = 1;
    }
  }
  fclose(trajfile);
}
