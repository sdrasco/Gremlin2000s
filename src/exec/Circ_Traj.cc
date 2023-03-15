//---------------------------------------------------------------------------
//
// $Id: Circ_Traj.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// This code will perform interpolations along a radiation reaction
// sequence to generate a trajectory through the parameter space of
// Circular Kerr geodesic orbits.
//
// Input parameters related to the trajectory's initial conditions and
// the files in which the radiation reaction quantities are stored,
// output the trajectory through (r, iota) space.
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
#include "CID.h"
#include "NRUtil.h"

#define DEG_TO_RAD (M_PI/180.0)

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  if(argc != 10) {
    cerr << "Arguments:  1. data file basename  2. imax" << endl;
    cerr << "            3. rmin  4. rmax  5. dr" << endl;
    cerr << "            6. rstart  7. iotastart  8. M  9. mu" << endl;
    exit(1);
  }
  //
  char basename[40], trajname[40];
  FILE *trajfile;
  //
  sprintf(basename, "%s", argv[1]);
  const int imax = atoi(argv[2]);
  const Real rmin = (Real)atof(argv[3]);
  const Real rmax = (Real)atof(argv[4]);
  const Real dr = (Real)atof(argv[5]);
  Real rstart = (Real)atof(argv[6]);
  Real iotastart = (Real)atof(argv[7]);
  const Real M = (Real)atof(argv[8]);
  const Real mu = (Real)atof(argv[9]);
  const Real eta = mu/M;
  //
  CID cid(basename, imax, rmin, rmax, dr, 0, 0);
  //
  // Start spiralling on in...
  //
  sprintf(trajname, "%s_r%2.1lf_i%d.traj", basename, rstart,
	  ((int)iotastart));
  if (!(trajfile = fopen(trajname, "w"))) {
    cerr << "Error opening " << trajname << endl;
    exit(3);
  }
  // bookkeeper
  int DONE_YET = 0;
  //
  // dt is chosen so that we get good resolution on the inspiral.
  //
  Real dt = eta/(10.*cid.Omega_max);
  //
  // Initial conditions
  //
  Real t = 0., N_ph_cyc = 0., N_th_cyc = 0.;
  Real rad = rstart, cosincl = cos(iotastart*DEG_TO_RAD);
  cid.Get_lso(cosincl);
  //
  // If we don't have at least two grid points, the splines are unhappy.
  //
  if (cid.r_lso > cid.r_arr[cid.jmax - 1])
    Die("Not enough data for this starting inclination!");
  while (!DONE_YET) {
    cid.Get_coordsdot_omegas(rad, cosincl);
    //
    //  1. (BL time)/M
    //  2. (radius)/M
    //  3. inclination angle
    //  4. Nphi
    //  5. Ntheta
    //  6. omega_phi
    //  7. omega_theta
    //
    fprintf(trajfile, "%e %lf %lf %lf %lf %lf %lf\n",
	    t, rad, acos(cosincl)/DEG_TO_RAD,
	    N_ph_cyc, N_th_cyc, cid.Omega_ph, cid.Omega_th);
    //
    // Ending condition: If dt gets way tiny, kill the loop.  We
    // decrease dt if the next step will go inside the minimum
    // radius or cross the LSO.
    //
    if (dt > 1.e-7) {
      Real dr_tmp = cid.rdot*dt;
      //
      // First, check whether we're likely to cross rmin.
      //
      while (rad + dr_tmp < cid.rmin && dt > 1.e-8) {
	dt /= 2.;
	dr_tmp = cid.rdot*dt;
      }
      //
      // Next, check whether we might cross the LSO:
      //
      Real dcosiota_tmp = cid.cosiotadot*dt;
      cid.Get_lso(cosincl);
      while (rad + dr_tmp < cid.r_lso && dt > 1.e-8) {
	dt /= 2.;
	dcosiota_tmp = cid.cosiotadot*dt;
	cid.Get_lso(cosincl + dcosiota_tmp);
	dr_tmp = cid.rdot*dt;
      }
      t += dt/eta;
      N_ph_cyc += cid.Omega_ph*dt/(2.*eta*M_PI);
      N_th_cyc += cid.Omega_th*dt/(2.*eta*M_PI);
      rad += dr_tmp;
      cosincl += dcosiota_tmp;
    } else {
      DONE_YET = 1;
    }
  }
  fclose(trajfile);
}
