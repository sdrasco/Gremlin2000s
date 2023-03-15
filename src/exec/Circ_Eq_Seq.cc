//---------------------------------------------------------------------------
//
// $Id: Circ_Eq_Seq.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
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
#include "CETD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  char basename[40], outname[50];
  if(argc != 9) {
    cerr << "Arguments: 1. r_start  2. r_stop  3. a" << endl;
    cerr << "           4. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           5. epsilon  6. lmax  7. dr" << endl;
    cerr << "           8. output file basename" << endl;
    exit(1);
  }
  const Real r_start = (Real)atof(argv[1]);
  const Real r_stop = (Real)atof(argv[2]);
  const Real a = (Real)atof(argv[3]);
  const int proret = atoi(argv[4]);
  Real r_isco;
  if (proret)
    r_isco = Kerr::isco_pro(a);
  else
    r_isco = Kerr::isco_ret(a);
  const Real r_min = Max(r_stop, r_isco);
  const Real EPSILON  = (Real)atof(argv[5]);
  const int lmax = atoi(argv[6]);
  const Real dr = (Real)atof(argv[7]);
  sprintf(basename, "%s", argv[8]);
  //
  Real r;
  for (r = r_start; r >= r_min; r -= dr) {
    if (dr < 0.1)
      sprintf(outname, "%s_r%3.2lf.bin", basename, r);
    else
      sprintf(outname, "%s_r%2.1lf.bin", basename, r);
    //
    CETD cetd(proret, r, a, EPSILON, outname);
    //
    cetd.Driver(lmax);
  }
}
