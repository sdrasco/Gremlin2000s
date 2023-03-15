//---------------------------------------------------------------------------
//
// $Id: Circ.cc,v 1.2 2004/09/04 19:28:21 sdrasco Exp $
//
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include "Globals.h"
#include "CKG.h"
#include "CTD.h"

int main(int argc, char **argv)
{
  ios::sync_with_stdio();
  char outname[50];
  //
  if(argc != 6) {
    cerr << "Arguments: 1. r  2. a  3. cosiota" << endl;
    cerr << "           4. eps  5. output file basename" << endl;
    exit(1);
  }
  const Real r = (Real)atof(argv[1]);
  const Real a = (Real)atof(argv[2]);
  const Real mu = (Real)atof(argv[3]);
  const Real eps = (Real)atof(argv[4]);
  sprintf(outname, "%s.bin", argv[5]);

  cout.precision(16);
  cout << "% r = " << r << endl
       << "% a = " << a << endl
       << "% cosiota = " << mu << endl
       << "% eps = " << eps << endl
       << "% below is: [l m k omega EH EInf LH LInf]" << endl;

  CTD ctd(r, mu, USING_cosiota, a, eps, outname);
  ctd.Driver(eps, eps);
}
