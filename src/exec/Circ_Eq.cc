//---------------------------------------------------------------------------
//
// $Id: Circ_Eq.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include "Globals.h"
#include "CETD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  char outname[50];
  //
  if(argc != 6) {
    cerr << "Arguments: 1. r  2. a" << endl;
    cerr << "           3. Prograde (1) or Retrograde (0)" << endl;
    cerr << "           4. eps  5. output file basename" << endl;
    exit(1);
  }

  const Real r = (Real)atof(argv[1]);
  const Real a = (Real)atof(argv[2]);
  const int proret = atoi(argv[3]);
  const Real eps = (Real)atof(argv[4]);
  sprintf(outname, "%s.bin", argv[5]);

  CETD cetd(proret, r, a, eps, outname);
  cetd.Driver(eps);
}
