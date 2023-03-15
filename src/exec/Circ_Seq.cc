//---------------------------------------------------------------------------
//
// $Id: Circ_Seq.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include "Globals.h"
#include "NRUtil.h"
#include "CKG.h"
#include "CTD.h"
#include "CETD.h"
#include "LB.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  char outname[30];
  //
  if(argc != 8) {
    cerr << "Arguments: 1. a  2. rmax  3. rmin  4. dr  5. N_i" << endl;
    cerr << endl;
    cerr << "   'i' is an integer index: i = 1 is the lowest energy" << endl;
    cerr << "    orbit at a given radius (prograde, equatorial);" << endl;
    cerr << "    i = N_i is the highest energy orbit (retrograde" << endl;
    cerr << "    equatorial, or orbit on separatrix)." << endl;
    cerr << endl;
    cerr << "           6. EPSILON  7. output base name" << endl;
    exit(1);
  }
  const Real a = (Real)atof(argv[1]);
  const Real rmax = (Real)atof(argv[2]);
  Real rmin = (Real)atof(argv[3]);
  const Real isco_pro = Kerr::isco_pro(a);
  if (rmin < isco_pro) rmin = isco_pro;
  const Real dr = (Real)atof(argv[4]);
  const int N_i = atoi(argv[5]);
  const Real EPSILON = (Real)atof(argv[6]);
  //
  int EQUATORIAL, RET_EXISTS, CURRENT_PRORET;
  int i;
  Real r;
  for (r = rmax; r >= rmin; r -= dr) {
    const Real Emb = Kerr::Eeqpro(r, a);
    const Real Lzmb = Kerr::Lzeqpro(r, a);
    Real Emax, Eninety_pro, Eninety_ret;
    LB lb(r, a);
    //
    // Do retrograde orbits exist at this radius?
    if (lb.cosiota < 0.) {
      RET_EXISTS = 1;
      //
      // Characteristics of the slightly prograde orbit
      CKG ckg_pro(1, 0, 1e-4, USING_cosiota, r, a);
      Eninety_pro = ckg_pro.E;
      Emax = Eninety_pro;
      //
      // Characteristics of the slightly retrograde orbit
      CKG ckg_ret(1, 0,-1e-4, USING_cosiota, r, a);
      Eninety_ret = ckg_ret.E;
    } else {
      RET_EXISTS = 0;
      Emax = lb.E;
    }
    //
    // Do the prograde sequence
    const Real deltaE = (Emax - Emb)/((Real)(N_i - 1));
    CURRENT_PRORET = 1;
    for (i = 1; i <= N_i; i++) {
      sprintf(outname, "%s_pro_r%2.1lf_i%d.bin", argv[7], r, i);
      Real E = Emb + ((Real)(i - 1))*deltaE;
      //
      // When E = Elb and the orbit is nonequatorial, the value of
      // Qdot looks rather wacky.  This isn't surprising since you are
      // right on the plunge orbit.  Move slightly off.
      //
      if (fabs(E - lb.E) < 1.e-10 && !EQUATORIAL)
	E -= deltaE*0.005;
      //
      // Is this an equatorial orbit?
      if (i == 1) EQUATORIAL = 1;
      else EQUATORIAL = 0;
      //
      if (EQUATORIAL) {
	CETD cetd(CURRENT_PRORET, r, a, Min(EPSILON, 1e-6), outname);
	//
	cetd.Driver(10. * EPSILON);
      } else {
	CTD ctd(r, E, USING_E, a, Min(EPSILON, 1e-6), outname);
	//
	ctd.Driver(10. * EPSILON, 10. * EPSILON);
      } // close 'else'
      cout << "DONE PROGRADE r = " << r << " " << "i = " << i << endl;
    } // close the i-loop
    //
    // Do the retrograde sequence, if necessary
    if (RET_EXISTS) {
      Emax = lb.E;
      const Real deltaE = (Emax - Eninety_ret)/((Real)(N_i - 1));
      CURRENT_PRORET = 0;
      for (i = 1; i <= N_i; i++) {
	sprintf(outname, "%s_ret_r%2.1lf_i%d.bin", argv[7], r, i);
	Real E = Eninety_ret + ((Real)(i - 1))*deltaE;
	//
	// When E = Elb and the orbit is nonequatorial, the value of
	// Qdot looks rather wacky.  This isn't surprising since you are
	// right on the plunge orbit.  Move slightly off.
	//
	if (fabs(E - lb.E) < 1.e-10 && !EQUATORIAL)
	  E -= deltaE*0.005;
	//
	// Is this an equatorial orbit?
	if (i == N_i && lb.Q == 0) EQUATORIAL = 1;
	else EQUATORIAL = 0;
	//
	if (EQUATORIAL) {
	  CETD cetd(CURRENT_PRORET, r, a, Min(EPSILON, 1e-6), outname);
	  //
	  cetd.Driver(10. * EPSILON);
	} else {
	  CTD ctd(r, E, USING_E, a, Min(EPSILON, 1e-6), outname);
	  //
	  ctd.Driver(10. * EPSILON, 10. * EPSILON);
	} // close 'else'
	cout << "DONE RETROGRADE r = " << r << " " << "i = " << i << endl;
      } // close the i-loop
    } // close if
  } // close the r-loop
}
