//---------------------------------------------------------------------------
//
// $Id: CEID.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular, equatorial inspiral data.
//
#include <math.h>
#include <stdio.h>
#include "Globals.h"
#include "CEID.h"
#include "SWSH.h"
#include "RRGW.h"
#include "NRUtil.h"

CEID::CEID(char *basename, const Real arrmin, const Real rmax, const Real dr,
	   const int ellmax) : rmin(arrmin), lmax(ellmax)
{
  char inname[50];
  //
  // This assignment of jmax is a bit whacky, but it works well.
  //
  jmax = (int)ceil((rmax - rmin)/dr + 0.0001);
  //
  // Allocate memory
  //
  Allocate();
  //
  // Read in data.  Note reversal of i index order --- want index to
  // increase in the same direction as cos(iota).
  //
  Omega_max = 0.;
  for (j = 1; j <= jmax; j++) {
    Real radius = rmin + (j - 1)*dr;
    if (dr < 0.1)
      sprintf(inname, "%s_r%3.2lf.bin", basename, radius);
    else
      sprintf(inname, "%s_r%2.1lf.bin", basename, radius);
    ReadIn(inname);
  }
  //
  // Make radial splines
  //
  Make_radial_splines();
}

void CEID::Allocate()
{
  r_arr = Realvector(1, jmax);
  rstar_arr = Realvector(1, jmax);
  Omega_arr = Realvector(1, jmax);
  rdot_arr = Realvector(1, jmax);
  Edot_arr = Realvector(1, jmax);
  Lzdot_arr = Realvector(1, jmax);
  //
  Omega_pp = Realvector(1, jmax);
  rdot_pp = Realvector(1, jmax);
  Edot_pp = Realvector(1, jmax);
  Lzdot_pp = Realvector(1, jmax);
  //
  if (lmax > 0) {
    Z_re = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    Z_re_pp = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    Z_im = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    Z_im_pp = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    Spheroid = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    Spheroid_pp = Real3tensor(2, lmax, -lmax, lmax, 1, jmax);
    //
    // Just to be on the safe side, these guys should all be
    // initialized to zero ... don't trust the compiler to do it
    // properly for me ...
    //
    for(j = 1; j <= jmax; j++) {
      for (l = 2; l <= lmax; l++) {
	for (m = -l; m <= l; m++) {
	  Z_re[l][m][j] = 0.;
	  Z_im[l][m][j] = 0.;
	  Spheroid[l][m][j] = 0.;
	}
      }
    }
    //
    Pplus_r = Realvector(1, jmax);
    Pplus_r_pp = Realvector(1, jmax);
    Pplus_i = Realvector(1, jmax);
    Pplus_i_pp = Realvector(1, jmax);
    Pminus_r = Realvector(1, jmax);
    Pminus_r_pp = Realvector(1, jmax);
    Pminus_i = Realvector(1, jmax);
    Pminus_i_pp = Realvector(1, jmax);
    Pz_r = Realvector(1, jmax);
    Pz_r_pp = Realvector(1, jmax);
    Pz_i = Realvector(1, jmax);
    Pz_i_pp = Realvector(1, jmax);
  }
}

void CEID::ReadIn(char *inname)
{
  FILE *infile;
  DataHolder indata;
  if (!(infile = fopen(inname, "r"))) {
    cerr << "File " << inname << " does not exist!!" << endl;
    exit(3);
  }
  //
  Real EdotInf = 0., EdotH = 0.;
  Real LzdotInf = 0., LzdotH = 0.;
  Real rdotInf = 0., rdotH = 0.;
  //
  while(fread(&indata, sizeof(DataHolder), 1, infile)) {
    a = indata.a;
    l = indata.l;
    m = indata.m;
    //
    // Factor of 2 because data only includes positive m, but we've
    // got symmetry.  Note that Edot and Lzdot are the fluxes
    // carried by the radiation, NOT the change in the particle's
    // energy and angular momentum.
    //
    EdotInf += 2.*indata.EdotInf;
    LzdotInf += 2.*indata.LzdotInf;
    rdotInf += 2.*indata.rdotInf;
    EdotH += 2.*indata.EdotH;
    LzdotH += 2.*indata.LzdotH;
    rdotH += 2.*indata.rdotH;
    //
    // Only store frequency for l == 2, m == 1
    //
    if (l == 2 && m == 1)
      Omega_arr[j] = indata.w;
    if (fabs(Omega_arr[j]) > Omega_max)
      Omega_max = fabs(Omega_arr[j]);
    //
    // ZH_{l, -m} = (-1)^l conj(ZH_{l,m}):
    //
    if (l <= lmax) {
      Z_re[l][m][j] = indata.ZH.real();
      Z_im[l][m][j] = indata.ZH.imag();
      Real sign = (l%2) ? -1. : 1.;
      Z_re[l][-m][j] = sign*Z_re[l][m][j];
      Z_im[l][-m][j] = -sign*Z_im[l][m][j];
    }
  } // while loop
  r_arr[j] = indata.r;
  rstar_arr[j] = Kerr::rstar(indata.r, a);
  fclose(infile);
  //
  Edot_arr[j] = EdotInf + EdotH;
  Lzdot_arr[j] = LzdotInf + LzdotH;
  //
  // Use r_isco to clear out divergence in rdot towards the lso
  //
  if (indata.Lz > 0.) { // prograde
    r_isco = Kerr::isco_pro(a);
  } else { // retrograde
    r_isco = Kerr::isco_ret(a);
  }
  rdot_arr[j] = (rdotInf + rdotH)*(r_arr[j] - r_isco);
}

void CEID::Make_radial_splines()
{
  spline(rstar_arr, Omega_arr, jmax, 1.e30, 1.e30, Omega_pp);
  spline(rstar_arr, rdot_arr, jmax, 1.e30, 1.e30, rdot_pp);
  spline(rstar_arr, Edot_arr, jmax, 1.e30, 1.e30, Edot_pp);
  spline(rstar_arr, Lzdot_arr, jmax, 1.e30, 1.e30, Lzdot_pp);
  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      spline(rstar_arr, Z_re[l][m], jmax, 1.e30, 1.e30, Z_re_pp[l][m]);
      spline(rstar_arr, Z_im[l][m], jmax, 1.e30, 1.e30, Z_im_pp[l][m]);
    }
  }
}

CEID::~CEID()
{
  free_Realvector(Lzdot_pp, 1, jmax);
  free_Realvector(Edot_pp, 1, jmax);
  free_Realvector(rdot_pp, 1, jmax);
  free_Realvector(Omega_pp, 1, jmax);
  //
  free_Realvector(Lzdot_arr, 1, jmax);
  free_Realvector(Edot_arr, 1, jmax);
  free_Realvector(rdot_arr, 1, jmax);
  free_Realvector(Omega_arr, 1, jmax);
  free_Realvector(rstar_arr, 1, jmax);
  free_Realvector(r_arr, 1, jmax);
  if (lmax > 0) {
    free_Real3tensor(Z_re, 2, lmax, -lmax, lmax, 1, jmax);
    free_Real3tensor(Z_re_pp, 2, lmax, -lmax, lmax, 1, jmax);
    free_Real3tensor(Z_im, 2, lmax, -lmax, lmax, 1, jmax);
    free_Real3tensor(Z_im_pp, 2, lmax, -lmax, lmax, 1, jmax);
    free_Real3tensor(Spheroid, 2, lmax, -lmax, lmax, 1, jmax);
    free_Real3tensor(Spheroid_pp, 2, lmax, -lmax, lmax, 1, jmax);
  }
}

void CEID::Spheroids(const Real costheta_view)
{
  if (lmax > 0) {
    for (j = 1; j <= jmax; j++) {
      for (l = 2; l <= lmax; l++) {
	for (m = -l; m <= l; m++) {
	  SWSH swsh(l, -2, m, a*m*Omega_arr[j]);
	  Spheroid[l][m][j] = swsh.spheroid(costheta_view);
	}
      }
    }
    for (l = 2; l <= lmax; l++)
      for (m = -l; m <= l; m++)
	spline(rstar_arr, Spheroid[l][m], jmax, 1.e30, 1.e30,
	       Spheroid_pp[l][m]);
  } else return;
}

void CEID::Pfluxes()
{
  if (lmax > 0) {
    for (j = 1; j <= jmax; j++) {
      Complex **Z; Z = Complexmatrix(2, lmax, -lmax, lmax);
      for (l = 2; l <= lmax; l++)
	for (m = -l; m <= l; m++)
	  Z[l][m] = Z_re[l][m][j] + II*Z_im[l][m][j];
      //
      Complex tmp_Pplus = Pplus(Z, Omega_arr[j])/(4.*M_PI);
      Pplus_r[j] = tmp_Pplus.real();
      Pplus_i[j] = tmp_Pplus.imag();
      Complex tmp_Pminus = Pminus(Z, Omega_arr[j])/(4.*M_PI);
      Pminus_r[j] = tmp_Pminus.real();
      Pminus_i[j] = tmp_Pminus.imag();
      Complex tmp_Pz = Pz(Z, Omega_arr[j])/(4.*M_PI);
      Pz_r[j] = tmp_Pz.real();
      Pz_i[j] = tmp_Pz.imag();
      //
      free_Complexmatrix(Z, 2, lmax, -lmax, lmax);
    }
    spline(rstar_arr, Pplus_r, jmax, 1.e30, 1.e30, Pplus_r_pp);
    spline(rstar_arr, Pplus_i, jmax, 1.e30, 1.e30, Pplus_i_pp);
    spline(rstar_arr, Pminus_r, jmax, 1.e30, 1.e30, Pminus_r_pp);
    spline(rstar_arr, Pminus_i, jmax, 1.e30, 1.e30, Pminus_i_pp);
    spline(rstar_arr, Pz_r, jmax, 1.e30, 1.e30, Pz_r_pp);
    spline(rstar_arr, Pz_i, jmax, 1.e30, 1.e30, Pz_i_pp);
  } else return;
}

void CEID::Get_rdot_omega(const Real r)
{
  const Real rs = Kerr::rstar(r, a);

  splint(rstar_arr, rdot_arr, rdot_pp, jmax, rs, &rdot);
  rdot /= (r - r_isco);
  splint(rstar_arr, Omega_arr, Omega_pp, jmax, rs, &Omega);
}

void CEID::Get_fluxes(const Real r)
{
  const Real rs = Kerr::rstar(r, a);

  splint(rstar_arr, Edot_arr, Edot_pp, jmax, rs, &Edot);
  splint(rstar_arr, Lzdot_arr, Lzdot_pp, jmax, rs, &Lzdot);
}

void CEID::Get_wave(const Real r, const Real phi, const Real phi_view)
{
  const Real rs = Kerr::rstar(r, a);

  hp = 0.; hc = 0.;
  Real Z_re_int, Z_im_int, Spheroid_int, dhc, dhp;
  Complex Zhere;
  RRGW rrgw;
  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      if (m != 0) {
	splint(rstar_arr, Z_re[l][m], Z_re_pp[l][m], jmax, rs, &Z_re_int);
	splint(rstar_arr, Z_im[l][m], Z_im_pp[l][m], jmax, rs, &Z_im_int);
	splint(rstar_arr, Spheroid[l][m], Spheroid_pp[l][m], jmax, rs,
	       &Spheroid_int);
	Zhere = Z_re_int + II*Z_im_int;
	rrgw.Wave(m, 0, phi/(2.*M_PI), 0., phi_view, Spheroid_int,
		  m*Omega, Zhere, dhp, dhc);
	hp += dhp; hc += dhc;
      }
    }
  }
}

// This code will not get your mother.  Sorry.
void CEID::Get_mom(const Real r, const Real phi)
{
  const Real rs = Kerr::rstar(r, a);

  Real Pplus_r_int, Pplus_i_int;
  Real Pminus_r_int, Pminus_i_int;
  Real Pz_r_int, Pz_i_int;
  //
  splint(rstar_arr, Pplus_r, Pplus_r_pp, jmax, rs, &Pplus_r_int);
  splint(rstar_arr, Pplus_i, Pplus_i_pp, jmax, rs, &Pplus_i_int);
  splint(rstar_arr, Pminus_r, Pminus_r_pp, jmax, rs, &Pminus_r_int);
  splint(rstar_arr, Pminus_i, Pminus_i_pp, jmax, rs, &Pminus_i_int);
  splint(rstar_arr, Pz_r, Pz_r_pp, jmax, rs, &Pz_r_int);
  splint(rstar_arr, Pz_i, Pz_i_pp, jmax, rs, &Pz_i_int);
  //
  Complex plus = exp(II*phi)*(Pplus_r_int + II*Pplus_i_int);
  Complex minus = exp(-II*phi)*(Pminus_r_int + II*Pminus_i_int);
  //
  // The next three combinations should all actually be _real_
  // numbers.  We leave them like this for now to make sure nothing is
  // acting screwy.  Note the sign --- we want system's recoil,
  // opposite of the momentum carried by the radiation.
  //
  pdot_x = -0.5*(plus + minus);
  pdot_y = II*0.5*(plus - minus);
  pdot_z = -(Pz_r_int + II*Pz_i_int);
}

Complex CEID::Pplus(Complex **Z, const Real omega)
{
  int l, lp, m, mp;
  Complex sum = 0.;
  for (l = 2; l <= lmax; l++) {
    for (lp = 2; lp <= lmax; lp++) {
      for (m = -l; m <= l; m++) {
	mp = m + 1;
	//
	// Note the values we exclude.  This is because
	// ZH_{lm} = 0 for m = 0
	// |m + 1| <= lp is required
	//
 	if (m != 0 && mp != 0 && abs(mp) <= lp) {
	  Complex Zed_l = Z[l][m];
	  Complex Zed_lp = Z[lp][m + 1];
	  
	  sum += Zed_l*conj(Zed_lp)*Il_lp_m_mp(l, lp, m, m + 1, omega)/
	    (m*(m + 1)*omega*omega);
	}
      }
    }
  }
  return(sum);
}

Complex CEID::Pminus(Complex **Z, const Real omega)
{
  int l, lp, m, mp;
  Complex sum = 0.;
  for (l = 2; l <= lmax; l++) {
    for (lp = 2; lp <= lmax; lp++) {
      for (m = -l; m <= l; m++) {
	mp = m - 1;
	//
	// Note the values we exclude.  This is because
	// ZH_{lm} = 0 for m = 0
	// |m - 1| <= lp is required
	//
 	if (m != 0 && mp != 0 && abs(mp) <= lp) {
	  Complex Zed_l = Z[l][m];
	  Complex Zed_lp = Z[lp][m - 1];

	  sum += Zed_l*conj(Zed_lp)*Il_lp_m_mp(l, lp, m, m - 1, omega)/
	    (m*(m - 1)*omega*omega);
 	}
      }
    }
  }
  return(sum);
}

Complex CEID::Pz(Complex **Z, const Real omega)
{
  int l, lp, m;
  Complex sum = 0;
  for (l = 2; l <= lmax; l++) {
    for (lp = 2; lp <= lmax; lp++) {
      for (m = -l; m <= l; m++) {
	//
	// Note the values we exclude.  This is because
	// ZH_{lm} = 0 for m = 0
	// |m| <= lp is required
	//
 	if (m != 0 && abs(m) <= lp) {
	  Complex Zed_l = Z[l][m];
	  Complex Zed_lp = Z[lp][m];

	  sum += Zed_l*conj(Zed_lp)*Jl_lp_m(l, lp, m, omega)/
	    (m*m*omega*omega);
 	}
      }
    }
  }
  if (abs(sum) < 1.e-15) sum = 0.;
  return(sum);
}

Real CEID::Il_lp_m_mp(const int l, const int lp, const int m, const int mp,
		      const Real omega)
{
  const Real aw = m*a*omega;
  const Real awp = mp*a*omega;
  //
  SWSH lm(l, -2, m, aw);
  SWSH lpmp(lp, -2, mp, awp);
  Clebsch cg;
  Real sum = 0.;
  int i, ip;
  int k, kp;
  //
  // Note some slightly complicated shenanigans required because the
  // way that the coefficients are stored.  b[1] corresponds to the
  // entry for lmin, b[2] corresponds to lmin + 1, etc.
  //
  for (i = 1; i <= lm.N; i++) {
    k = i - 1 + lm.lmin;
    //
    // First term
    //
    kp = k - 1;
    // Find corresponding index.
    ip = kp + 1 - lpmp.lmin;
    if (ip > 0 && ip <= lpmp.N) {
      sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
    }
    //
    // Second term
    //
    kp = k;
    ip = kp + 1 - lpmp.lmin;
    if (ip > 0 && ip <= lpmp.N) {
      sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
    }
    //
    // Third term
    //
    kp = k + 1;
    ip = kp + 1 - lpmp.lmin;
    if (ip > 0 && ip <= lpmp.N) {
      sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
    }
  }
  return(sum);
}

Real CEID::Jl_lp_m(const int l, const int lp, const int m, const Real omega)
{
  const Real aw = m*a*omega;
  //
  SWSH lm(l, -2, m, aw);
  SWSH lpm(lp, -2, m, aw);
  Clebsch cg;
  Real sum = 0.;
  int i, ip;
  int k, kp;
  //
  // Note some slightly complicated shenanigans required because the
  // way that the coefficients are stored.  b[1] corresponds to the
  // entry for lmin, b[2] corresponds to lmin + 1, etc.
  //
  for (i = 1; i <= lm.N; i++) {
    k = i - 1 + lm.lmin;
    //
    // First term
    //
    kp = k - 1;
    // Find corresponding index.
    ip = kp + 1 - lpm.lmin;
    if (ip > 0 && ip <= lpm.N) {
      sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
    }
    //
    // Second term
    //
    kp = k;
    ip = kp + 1 - lpm.lmin;
    if (ip > 0 && ip <= lpm.N) {
      sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
    }
    //
    // Third term
    //
    kp = k + 1;
    ip = kp + 1 - lpm.lmin;
    if (ip > 0 && ip <= lpm.N) {
      sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
    }
  }
  return(sum);
}
