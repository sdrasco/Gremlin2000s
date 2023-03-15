//---------------------------------------------------------------------------
//
// $Id: GW.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
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
#include "SWSH.h"
#include "RRGW.h"
#define SRATE 12000.0      // Hz
#define MSUN 4.92e-6       // seconds
#define HOUR 3600.
#define CLIGHT 3.00e10     // cm/sec
#define MEGAPARSEC 3.09e24 // cm

int main(int argc, char **argv)
{
  ios::sync_with_stdio();

  int i, l, m, k, harm_count = 0;
  const int lmax = 30, krange = 40;
  Real EdotNewt(const Real r);
  Real E, Lz, Q;
  Real r, a, cosiota, iotaprime, t, t_sec, t_hr, t_day, tmincyc, dt;
  Real dhp, dhc;
  Real hp_phys, hc_phys;
  Real hp_geom, hc_geom;
  Real wmin = 1.e30, wmax = 0.;
  Real **w, **p;
  Real *w2;
  //
  // Word of caution: Edot and Lzdot in the datafile are the
  // energy and angular momentum rate of change in the radiation.
  // The change in the particle's energy and angular momentum is
  // Edot_particle = -Edot, Lzdot_particle = -Lzdot.
  //
  // Qdot and rdot in the datafile are the rates of change for
  // the particle.
  //
  Real ***EdotInf, ***LzdotInf, ***QdotInf, ***rdotInf;
  Real **EdotInf_mk, **LzdotInf_mk, **QdotInf_mk, **rdotInf_mk;
  Real ***EdotH, ***LzdotH, ***QdotH, ***rdotH;
  Real **EdotH_mk, **LzdotH_mk, **QdotH_mk, **rdotH_mk;
  Real *EdotH_w, *EdotInf_w;
  Real EdotInf_tot = 0., LzdotInf_tot = 0., QdotInf_tot = 0.,
    rdotInf_tot = 0.;
  Real EdotH_tot = 0., LzdotH_tot = 0., QdotH_tot = 0., rdotH_tot = 0.;
  Real Edot_tot = 0., Lzdot_tot = 0., Qdot = 0., rdot = 0.;
  Real Edot_particle, Lzdot_particle;
  Real cosiotadot, iotadot;
  Complex ***ZH, ***ZInf;
  Real ***Spheroid;
  RRGW rrgw;
  char inname[40], specbasename[40], specname[40], sumspecname[40],
    wavename[40];
  DataHolder indata;
  FILE *infile, *specfile, *sumspecfile, *wavefile;

  if(argc != 5 && argc != 6 && argc != 8) {
    cerr << "Arguments: 1. Input file basename.  2. cos(theta)." << endl;
    cerr << "           3. phi (degrees).                      " << endl;
    cerr << "           4. Number of wmin cycles to go through." << endl;
    cerr << "(optional) 5. Black hole mass (in solar masses).  " << endl;
    cerr << "              If not given, it is set to 1 Msun.  " << endl;
    cerr << "(optional) 6. Particle mass (in solar masses).    " << endl;
    cerr << "(optional) 7. Distance to source (in Mpc).        " << endl;
    cerr << "              Args 5 & 6 must be given together.  " << endl;
    exit(1);
  }
  sprintf(inname, "%s.bin", argv[1]);
  sprintf(specbasename, "%s", argv[1]);
  sprintf(wavename, "%s.wave", argv[1]);
  const Real ct = (Real)atof(argv[2]);
  const Real phi = M_PI*(Real)atof(argv[3])/180.;
  const Real Nmincyc = atof(argv[4]);
  Real BHMass;
  if (argc > 5)
    BHMass = (Real)atof(argv[5]);
  else
    BHMass = 1.;
  Real PMass, SourceR, Amplitude;
  if (argc == 8) {
    PMass = (Real)atof(argv[6]);
    SourceR = (Real)atof(argv[7]);

    Amplitude = PMass*MSUN/(SourceR*MEGAPARSEC/CLIGHT);
  } else
    Amplitude = 1;
  //
  // Allocate memory.
  //
  w = Realmatrix(-lmax, lmax, -krange, krange);
  p = Realmatrix(-lmax, lmax, -krange, krange);
  EdotInf = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  LzdotInf = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  QdotInf = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  rdotInf = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  EdotH = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  LzdotH = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  QdotH = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  rdotH = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  EdotInf_mk = Realmatrix(-lmax, lmax, -krange, krange);
  LzdotInf_mk = Realmatrix(-lmax, lmax, -krange, krange);
  QdotInf_mk = Realmatrix(-lmax, lmax, -krange, krange);
  rdotInf_mk = Realmatrix(-lmax, lmax, -krange, krange);
  EdotH_mk = Realmatrix(-lmax, lmax, -krange, krange);
  LzdotH_mk = Realmatrix(-lmax, lmax, -krange, krange);
  QdotH_mk = Realmatrix(-lmax, lmax, -krange, krange);
  rdotH_mk = Realmatrix(-lmax, lmax, -krange, krange);
  ZH = Complex3tensor(2, lmax, -lmax, lmax, -krange, krange);
  ZInf = Complex3tensor(2, lmax, -lmax, lmax, -krange, krange);
  Spheroid = Real3tensor(2, lmax, -lmax, lmax, -krange, krange);
  //
  // Pre-load these guys with zeros; this way there is no
  // blank memory for the multipoles that don't get computed.
  //
  for(l = 2; l <= lmax; l++) {
    for(m = -l; m <= l; m++) {
      for(k = -krange; k <= krange; k++) {
	w[m][k] = 0.;
	EdotInf[l][m][k] = 0.;
	LzdotInf[l][m][k] = 0.;
	QdotInf[l][m][k] = 0.;
	rdotInf[l][m][k] = 0.;
	EdotH[l][m][k] = 0.;
	LzdotH[l][m][k] = 0.;
	QdotH[l][m][k] = 0.;
	rdotH[l][m][k] = 0.;
	ZH[l][m][k] = Complex(0., 0.);
	ZInf[l][m][k] = Complex(0., 0.);
      }
    }
  }
  //
  // Read in.
  //
  if (!(infile = fopen(inname, "r"))) {
    cerr << inname << " does not exist!" << endl;
    exit(0);
  }
  while(!feof(infile)) {
    fread(&indata, sizeof(DataHolder), 1, infile);
    l = indata.l;
    m = indata.m;
    k = indata.k;
    w[m][k] = indata.w;
    if(fabs(indata.w) > wmax)
      wmax = fabs(indata.w);
    if(fabs(indata.w) < wmin && fabs(indata.w) > 1.e-7)
      wmin = fabs(indata.w);
    p[m][k] = indata.p;
    ZH[l][m][k] = indata.ZH;
    ZInf[l][m][k] = indata.ZInf;
    EdotInf[l][m][k] = indata.EdotInf;
    LzdotInf[l][m][k] = indata.LzdotInf;
    QdotInf[l][m][k] = indata.QdotInf;
    rdotInf[l][m][k] = indata.rdotInf;
    EdotH[l][m][k] = indata.EdotH;
    LzdotH[l][m][k] = indata.LzdotH;
    QdotH[l][m][k] = indata.QdotH;
    rdotH[l][m][k] = indata.rdotH;
    if(EdotInf[l][m][k] != 0.0) harm_count++;
  }
  E = indata.E;
  a = indata.a;
  Lz = indata.Lz;
  Q = indata.Q;
  r = indata.r;
  int Nw = (int)(wmax/wmin) + 1;
  w2 = Realvector(1, Nw);
  EdotInf_w = Realvector(1, Nw);
  EdotH_w = Realvector(1, Nw);
  const Real Delta = Kerr::Delta(r, a);
  cosiota = Lz/sqrt(Lz*Lz + Q);
  iotaprime = atan2(sqrt(Q), Lz + a*(2.*E*r - a*Lz)/Delta);
  fprintf(stderr,"%d harmonics computed for this orbit.\n",harm_count);
  //
  // Use symmetry to fill in blank spots.
  //
  for(l = 2; l <= lmax; l++) {
    for(m = -l; m <= l; m++) {
      for(k = -krange; k <= krange; k++) {
	if (w[m][k] == 0.) w[m][k] = -w[-m][-k];
	if (EdotInf[l][m][k] == 0.) EdotInf[l][m][k] = EdotInf[l][-m][-k];
	if (LzdotInf[l][m][k] == 0.) LzdotInf[l][m][k] = LzdotInf[l][-m][-k];
	if (EdotH[l][m][k] == 0.) EdotH[l][m][k] = EdotH[l][-m][-k];
	if (LzdotH[l][m][k] == 0.) LzdotH[l][m][k] = LzdotH[l][-m][-k];
	if (QdotInf[l][m][k] == 0.) QdotInf[l][m][k] = QdotInf[l][-m][-k];
	if (QdotH[l][m][k] == 0.) QdotH[l][m][k] = QdotH[l][-m][-k];
	if (rdotInf[l][m][k] == 0.) rdotInf[l][m][k] = rdotInf[l][-m][-k];
	if (rdotH[l][m][k] == 0.) rdotH[l][m][k] = rdotH[l][-m][-k];
	Real sign;
	if (k%2)
	  sign = -1.;
	else
	  sign = 1.;
	if (ZInf[l][m][k] == Complex(0., 0.))
	  ZInf[l][m][k] = sign*conj(ZInf[l][-m][-k]);
	if (ZH[l][m][k] == Complex(0., 0.))
	  ZH[l][m][k] = sign*conj(ZH[l][-m][-k]);
      }
    }
  }
  fclose(infile);
  //
  // Write the spectral info.
  //
  double leinf;
  for(l = 2; l <= lmax; l++) {
    sprintf(specname, "%s_l%d.spec", specbasename, l);
    specfile = fopen(specname, "w");
    //
    // Find minimum flux to infinity on file...
    //
    double Edotmin = 666.;
    for(m = -l; m <= l; m++) {
      for(k = -krange+10; k <= krange-10; k++) {
	if (EdotInf[l][m][k] < Edotmin && EdotInf[l][m][k] > 0.0)
	  Edotmin = EdotInf[l][m][k];
      }
    }
    for(m = -l; m <= l; m++) {
      for(k = -krange+10; k <= krange-10; k++) {
	if (EdotInf[l][m][k] != 0.0) leinf = log10(EdotInf[l][m][k]);
	else leinf = log10(Edotmin);
     	fprintf(specfile, "%d %d %e %e %e %e %e %e %e %e %e %e\n",
     		m, k, w[m][k],
     		EdotInf[l][m][k], EdotH[l][m][k], leinf,
     		LzdotInf[l][m][k], LzdotH[l][m][k],
     		QdotInf[l][m][k], QdotH[l][m][k],
     		rdotInf[l][m][k], rdotH[l][m][k]);
	EdotInf_tot += EdotInf[l][m][k];
	LzdotInf_tot += LzdotInf[l][m][k];
	QdotInf_tot += QdotInf[l][m][k];
	rdotInf_tot += rdotInf[l][m][k];
	EdotH_tot += EdotH[l][m][k];
	LzdotH_tot += LzdotH[l][m][k];
	QdotH_tot += QdotH[l][m][k];
	rdotH_tot += rdotH[l][m][k];
      }
    }
    fclose(specfile);
  }
  for (i = 1; i <= Nw; i++) {
    w2[i] = (i - 1)*wmax/Nw;
    EdotInf_w[i] = 0.0;
    EdotH_w[i] = 0.0;
  }
  for (l = 2; l <= lmax; l++ ) {
    for (m = -l; m <= l; m++) {
      for (k = -krange; k <= krange; k++) {
	for (i = 1; i <= Nw; i++) {
	  if (w[m][k] > w2[i] && w[m][k] < w2[i+1]) {
	    EdotInf_w[i] += 2.*EdotInf[l][m][k];
	    EdotH_w[i] += 2.*EdotH[l][m][k];
	  }
	}
      }
    }
  }
  sprintf(sumspecname, "%s.sumspec", specbasename);
  sumspecfile = fopen(sumspecname, "w");
  for (i = 1; i <= Nw; i++)
    fprintf(sumspecfile, "%e %e %e\n", w2[i], EdotH_w[i], EdotInf_w[i]);
  fclose(sumspecfile);
  Edot_tot = EdotInf_tot + EdotH_tot;
  Edot_particle = -Edot_tot;
  Lzdot_tot = LzdotInf_tot + LzdotH_tot;
  Lzdot_particle = -Lzdot_tot;
  Qdot = QdotInf_tot + QdotH_tot;
  rdot = rdotInf_tot + rdotH_tot;
  cosiotadot = (Lzdot_particle - 0.5*Lz*(2.*Lz*Lzdot_particle + Qdot)/
		(Lz*Lz + Q))/sqrt(Lz*Lz + Q);
  iotadot = -cosiotadot/sqrt(1. - cosiota*cosiota);
  //
  // Write the waveform.
  //
  // First, do all the spheroidal harmonics.  They need only be
  // calculated once and stored.
  //
  for(l = 2; l <= lmax; l++) {
    for(m = -l; m <= l; m++) {
      for(k = -krange; k <= krange; k++) {
	if(w[m][k] != 0.0) {
	  SWSH swsh(l, -2, m, a*w[m][k]);
	  Spheroid[l][m][k] = swsh.spheroid(ct);
	}
	else 
	  Spheroid[l][m][k] = 0.0;
      }
    }
  }
  //
  // Do the high resolution waveform, useful for plots.
  //
  wavefile = fopen(wavename, "w");
  tmincyc = 2.*M_PI/wmin;
  dt = 2.*M_PI/(20.*wmax);
  //
  // First, some useful header info.
  //
  fprintf(wavefile, "# Omega_phi = %e  Omega_theta = %e\n",
	  w[1][0], w[0][1]);
  fprintf(wavefile, "# E = %e  Lz = %e  Q = %e  r = %e\n",
	  E, Lz, Q, r);
  fprintf(wavefile, "# iota      = %e (radians) %e (degrees)\n",
	  acos(cosiota), 180.*acos(cosiota)/M_PI);
  fprintf(wavefile, "# iotaprime = %e (radians) %e (degrees)\n",
	  iotaprime, 180.*iotaprime/M_PI);
  fprintf(wavefile, "# EdotInf = %e  LzdotInf = %e\n",
	  EdotInf_tot, LzdotInf_tot);
  fprintf(wavefile, "# QdotInf = %e  rdotInf = %e\n",
	  QdotInf_tot, rdotInf_tot);
  fprintf(wavefile, "# EdotH = %e  LzdotH = %e\n",
	  EdotH_tot, LzdotH_tot);
  fprintf(wavefile, "# QdotH = %e  rdotH = %e\n",
	  QdotH_tot, rdotH_tot);
  fprintf(wavefile, "# The following quantities are for the particle.\n");
  fprintf(wavefile, "# Edot = %e  Lzdot = %e  Qdot = %e  rdot = %e\n",
	  Edot_particle, Lzdot_particle, Qdot, rdot);
  fprintf(wavefile, "# cosiotadot = %e\n", cosiotadot);
  fprintf(wavefile, "# iotadot = %e\n", iotadot);
  //
  // Then, the waveform itself.
  //
  for(t = 0.; t <= Nmincyc*tmincyc; t += dt) {
    t_sec = t * BHMass * MSUN;
    t_hr = t_sec / HOUR;
    t_day = t_hr / 24.;
    hp_phys = 0.; hc_phys = 0.;
    hp_geom = 0.; hc_geom = 0.;
    for(l = 2; l <= lmax; l++) {
      for(m = -l; m <= l; m++) {
	for(k = -krange; k <= krange; k++) {
	  if(w[m][k] != 0 && (m != 0 || k != 0)) {
	    rrgw.Wave(m, t, phi, Spheroid[l][m][k],
		      w[m][k], ZH[l][m][k], dhp, dhc);
	    hp_phys += Amplitude*dhp; hc_phys += Amplitude*dhc;
	    hp_geom += dhp; hc_geom += dhc;
	  }
	}
      }
    }
    fprintf(wavefile, "%e %e %e %e %e %e %e %e\n", t, t_sec, t_hr, t_day,
	    hp_phys, hc_phys, hp_geom, hc_geom);
  }
  fclose(wavefile);
  //
  // Free memory.
  //
  free_Realmatrix(w, -lmax, lmax, -krange, krange);
  free_Realmatrix(p, -lmax, lmax, -krange, krange);
  free_Real3tensor(EdotInf, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(LzdotInf, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(QdotInf, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(rdotInf, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(EdotH, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(LzdotH, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(QdotH, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(rdotH, 2, lmax, -lmax, lmax, -krange, krange);
  free_Complex3tensor(ZH, 2, lmax, -lmax, lmax, -krange, krange);
  free_Complex3tensor(ZInf, 2, lmax, -lmax, lmax, -krange, krange);
  free_Real3tensor(Spheroid, 2, lmax, -lmax, lmax, -krange, krange);
  free_Realmatrix(EdotInf_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(LzdotInf_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(QdotInf_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(rdotInf_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(EdotH_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(LzdotH_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(QdotH_mk, -lmax, lmax, -krange, krange);
  free_Realmatrix(rdotH_mk, -lmax, lmax, -krange, krange);
}

Real EdotNewt(const Real r)
{
  return(6.4/(r*r*r*r*r));
}
