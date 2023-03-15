//---------------------------------------------------------------------------
//
// $Id: CID.cc,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Methods for circular inspiral data.
//
#include <math.h>
#include <stdio.h>
#include "Globals.h"
#include "CID.h"
#include "SWSH.h"
#include "LB.h"
#include "RRGW.h"
#include "NRUtil.h"

CID::CID(char *basename, const int ayemax, const Real arrmin,
	 const Real rmax, const Real dr, const int ellmax, const int kaymax)
  : rmin(arrmin), imax(ayemax), lmax(ellmax), kmax(kaymax)
{
  char inname[50];
  //
  // This assignment of jmax is a bit whacky, but it works well.
  //
  jmax = (int)ceil((rmax - rmin)/dr + 0.0001);
  //
  // Allocate memory.
  //
  Allocate();
  //
  // Read in data
  //
  Omega_max = 0.;
  for (j = 1; j <= jmax; j++) {
    Real radius = rmin + (j - 1)*dr;
    for (i = 1; i <= imax; i++) {
      sprintf(inname, "%s_pro_r%2.1f_i%d.bin", basename, radius, i);
      ReadIn(inname, PROGRADE);
      sprintf(inname, "%s_ret_r%2.1f_i%d.bin", basename, radius, i);
      ReadIn(inname, RETROGRADE);
    } // close i-loop
  } // close j-loop
  //
  // Make data needed for LSO.
  //
  Make_lso();
  //
  // Make angular splines.
  //
  Make_angular_splines();
}

void CID::Allocate()
{
  ret_exists = ivector(1, jmax);
  //
  r_arr = Realvector(1, jmax);
  rstar_arr = Realvector(1, jmax);
  //
  cosiota_arr = Realmatrix(1, jmax, 1, 2*imax);
  Omega_ph_arr = Realmatrix(1, jmax, 1, 2*imax);
  Omega_th_arr = Realmatrix(1, jmax, 1, 2*imax);
  rdot_arr = Realmatrix(1, jmax, 1, 2*imax);
  cosiotadot_arr = Realmatrix(1, jmax, 1, 2*imax);
  Edot_arr = Realmatrix(1, jmax, 1, 2*imax);
  Lzdot_arr = Realmatrix(1, jmax, 1, 2*imax);
  Qdot_arr = Realmatrix(1, jmax, 1, 2*imax);
  //
  Omega_ph_ipp = Realmatrix(1, jmax, 1, 2*imax);
  Omega_th_ipp = Realmatrix(1, jmax, 1, 2*imax);
  rdot_ipp = Realmatrix(1, jmax, 1, 2*imax);
  cosiotadot_ipp = Realmatrix(1, jmax, 1, 2*imax);
  Edot_ipp = Realmatrix(1, jmax, 1, 2*imax);
  Lzdot_ipp = Realmatrix(1, jmax, 1, 2*imax);
  Qdot_ipp = Realmatrix(1, jmax, 1, 2*imax);
  //
  // Just to be on the safe side, these guys should all be initialized
  // to zero ... don't trust the compiler to do it properly for me ...
  //
  for (j = 1; j <= jmax; j++) {
    for (i = 1; i <= 2*imax; i++) {
      Omega_ph_arr[j][i] = 0.;
      Omega_th_arr[j][i] = 0.;
      rdot_arr[j][i] = 0.;
      cosiotadot_arr[j][i] = 0.;
      Edot_arr[j][i] = 0.;
      Lzdot_arr[j][i] = 0.;
      Qdot_arr[j][i] = 0.;
      //
      Omega_ph_ipp[j][i] = 0.;
      Omega_th_ipp[j][i] = 0.;
      rdot_ipp[j][i] = 0.;
      cosiotadot_ipp[j][i] = 0.;
      Edot_ipp[j][i] = 0.;
      Lzdot_ipp[j][i] = 0.;
      Qdot_ipp[j][i] = 0.;
    }
  }
  //
  r_lso_arr = Realvector(1, LSO_POINTS);
  cosiota_lso_arr = Realvector(1, LSO_POINTS);
  r_lso_ipp = Realvector(1, LSO_POINTS);
  //
  Omega_ph_r = Realvector(1, jmax);
  Omega_th_r = Realvector(1, jmax);
  rdot_r = Realvector(1, jmax);
  cosiotadot_r = Realvector(1, jmax);
  Edot_r = Realvector(1, jmax);
  Lzdot_r = Realvector(1, jmax);
  Qdot_r = Realvector(1, jmax);
  //
  Omega_ph_rpp = Realvector(1, jmax);
  Omega_th_rpp = Realvector(1, jmax);
  rdot_rpp = Realvector(1, jmax);
  cosiotadot_rpp = Realvector(1, jmax);
  Edot_rpp = Realvector(1, jmax);
  Lzdot_rpp = Realvector(1, jmax);
  Qdot_rpp = Realvector(1, jmax);
  //
  if (lmax > 0) {
    Z_re = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax, 1, 2*imax);
    Z_im = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax, 1, 2*imax);
    Spheroid = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
			   1, 2*imax);
    //
    Z_re_ipp = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
			   1, 2*imax);
    Z_im_ipp = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
			   1, 2*imax);
    Spheroid_ipp = Real5tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
			       1, 2*imax);
    //
    Z_re_r = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    Z_im_r = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    Spheroid_r = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    //
    Z_re_rpp = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    Z_im_rpp = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    Spheroid_rpp = Real4tensor(2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    //
    // Just to be on the safe side, these guys should all be initialized
    // to zero ... don't trust the compiler to do it properly for me ...
    //
    for (l = 2; l <= lmax; l++) {
      for (m = -lmax; m <= lmax; m++) {
	for (k = -kmax; k <= kmax; k++) {
	  for (j = 1; j <= jmax; j++) {
	    Z_re_r[l][m][k][j] = 0.;
	    Z_im_r[l][m][k][j] = 0.;
	    Spheroid_r[l][m][k][j] = 0.;
	    Z_re_rpp[l][m][k][j] = 0.;
	    Z_im_rpp[l][m][k][j] = 0.;
	    Spheroid_rpp[l][m][k][j] = 0.;
	    for (i = 1; i <= 2*imax; i++) {
	      Z_re[l][m][k][j][i] = 0.;
	      Z_im[l][m][k][j][i] = 0.;
	      Spheroid[l][m][k][j][i] = 0.;
	      Z_re_ipp[l][m][k][j][i] = 0.;
	      Z_im_ipp[l][m][k][j][i] = 0.;
	      Spheroid_ipp[l][m][k][j][i] = 0.;
	    } // j-loop
	  } // i-loop
	} // k-loop
      } // m-loop
    } // l-loop
  } // if (lmax > 0)
}

void CID::ReadIn(char *inname, const int proflag)
{
  FILE *infile;
  DataHolder indata;
  int offset;
  //
  if (proflag) {
    if (!(infile = fopen(inname, "r"))) {
      cerr << "File " << inname << " does not exist!!" << endl;
      exit(3);
    } else {
      //
      // When reversing, data will live in the arrays from imax + 1 to
      // 2 imax.  Set offset to 2 imax in preparation.
      //
      offset = 2*imax;
    }
  } else { // !proflag
    if (!(infile = fopen(inname, "r"))) {
      //
      // Retrograde doesn't exist here.
      //
      ret_exists[j] = 0;
      return;
    } else {
      ret_exists[j] = 1;
      //
      // When reversing, data will live in the arrays from 1 to imax.
      // Set offset to imax in preparation.
      //
      offset = imax;
    }
  }
  Real EdotInf = 0., EdotH = 0.;
  Real LzdotInf = 0., LzdotH = 0.;
  Real QdotInf = 0., QdotH = 0.;
  Real rdotInf = 0., rdotH = 0.;
  while (fread(&indata, sizeof(DataHolder), 1, infile)) {
    a = indata.a;
    l = indata.l;
    m = indata.m;
    k = indata.k;
    //
    // Factor of 2 because data only includes positive k, but we've
    // got symmetry.
    //
    EdotInf += 2.*indata.EdotInf;
    LzdotInf += 2.*indata.LzdotInf;
    QdotInf += 2.*indata.QdotInf;
    rdotInf += 2.*indata.rdotInf;
    EdotH += 2.*indata.EdotH;
    LzdotH += 2.*indata.LzdotH;
    QdotH += 2.*indata.QdotH;
    rdotH += 2.*indata.rdotH;
    //
    // Only store frequency for l == 2.
    //
    if (l == 2) {
      if (m == 0 && k == 1)
	Omega_th_arr[j][offset - i + 1] = indata.w;
      if (m == 1 && k == 0) {
	Omega_ph_arr[j][offset - i + 1] = indata.w;
	if (indata.Q == 0.) {
	  //
	  // Equator: compute omega_th by hand (rr code puts zero
	  // since it doesn't compute any k dependence there!).
	  //
	  CKG ckg(0, 1, indata.E, USING_E, indata.r, a);
	  Omega_th_arr[j][offset - i + 1] = ckg.wmk;
	}
      } // close 'if (m == 1 && k == 0)'
      if (fabs(Omega_ph_arr[j][offset - i + 1]) > Omega_max)
	Omega_max = fabs(Omega_ph_arr[j][offset - i + 1]);
    } // close 'if (l == 2)'
    //
    // ZH_{l,-m,-k} = (-1)^(l + k) conj(ZH_{l,m,k}):
    //
    if (l <= lmax && abs(k) <= kmax) {
      Z_re[l][m][k][j][offset - i + 1] = indata.ZH.real();
      Z_im[l][m][k][j][offset - i + 1] = indata.ZH.imag();
      Real sign = ((l + k)%2) ? -1. : 1.;
      Z_re[l][-m][-k][j][offset - i + 1] = sign*Z_re[l][m][k][j][offset - i
								 + 1];
      Z_im[l][-m][-k][j][offset - i + 1] = -sign*Z_im[l][m][k][j][offset - i
								  + 1];
    }
  } // close while loop
  Real Lz = indata.Lz;
  Real Q = indata.Q;
  cosiota_arr[j][offset - i + 1] = Lz/sqrt(Lz * Lz + Q);
  r_arr[j] = indata.r;
  rstar_arr[j] = Kerr::rstar(indata.r, a);
  //
  fclose(infile);
  //
  Edot_arr[j][offset - i + 1] = EdotInf + EdotH;
  Lzdot_arr[j][offset - i + 1] = LzdotInf + LzdotH;
  //
  // This Lzdot is only used to compute cosiotadot --- switch sign so
  // that it applies to the particle.
  //
  Real Lzdot = -(LzdotInf + LzdotH);
  Qdot_arr[j][offset - i + 1] = QdotInf + QdotH;
  Real Qdot = Qdot_arr[j][offset - i + 1];
  cosiotadot_arr[j][offset - i + 1] = 1./sqrt(Lz*Lz + Q)*
    (Lzdot - 0.5*Lz*(2.*Lz*Lzdot + Qdot)/(Lz*Lz + Q));
  if (fabs(cosiotadot_arr[j][offset - i + 1]) < 1.e-14)
    cosiotadot_arr[j][offset - i + 1] = 0.;
  rdot_arr[j][offset - i + 1] = rdotInf + rdotH;
}

void CID::Make_lso()
{
  Real isco_pro = Kerr::isco_pro(a);
  Real isco_ret = Kerr::isco_ret(a);
  for (i = 1; i <= LSO_POINTS; i++) {
    r_lso_arr[i] = isco_ret - (i - 1)*(isco_ret - isco_pro)/
      ((Real)(LSO_POINTS - 1));
    LB lb(r_lso_arr[i], a);
    cosiota_lso_arr[i] = lb.cosiota;
  }
  spline(cosiota_lso_arr, r_lso_arr, LSO_POINTS, 0., 0., r_lso_ipp);
}

void CID::Make_angular_splines()
{
  for (j = 1; j <= jmax; j++) {
    for (i = 1; i <= 2*imax; i++) {
      Get_lso(cosiota_arr[j][i]);
      rdot_arr[j][i] *= (r_arr[j] - r_lso);
    }
    if (ret_exists[j]) {
      Ang_Spline(cosiota_arr[j], rdot_arr[j], rdot_ipp[j], FULL_RANGE);
      Ang_Spline(cosiota_arr[j], cosiotadot_arr[j], cosiotadot_ipp[j],
		 FULL_RANGE);
      Ang_Spline(cosiota_arr[j], Edot_arr[j], Edot_ipp[j], FULL_RANGE);
      Ang_Spline(cosiota_arr[j], Lzdot_arr[j], Lzdot_ipp[j], FULL_RANGE);
      Ang_Spline(cosiota_arr[j], Qdot_arr[j], Qdot_ipp[j], FULL_RANGE);
      Ang_Spline(cosiota_arr[j], Omega_th_arr[j], Omega_th_ipp[j], FULL_RANGE);
      //
      // Break Omega_ph into subranges.
      //
      Ang_Spline(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j],
		 UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j],
		 LOWER_RANGE);
      for (l = 2; l <= lmax; l++) {
	for (m = -lmax; m <= lmax; m++) {
	  for (k = -kmax; k <= kmax; k++) {
	    if (m == 0) {
	      Ang_Spline(cosiota_arr[j], Z_re[l][m][k][j],
			 Z_re_ipp[l][m][k][j], FULL_RANGE);
	      Ang_Spline(cosiota_arr[j], Z_im[l][m][k][j],
			 Z_im_ipp[l][m][k][j], FULL_RANGE);
	    } else {
	      //
	      // Break into subranges
	      //
	      Ang_Spline(cosiota_arr[j], Z_re[l][m][k][j],
			 Z_re_ipp[l][m][k][j], UPPER_RANGE);
	      Ang_Spline(cosiota_arr[j], Z_im[l][m][k][j],
			 Z_im_ipp[l][m][k][j], UPPER_RANGE);
	      Ang_Spline(cosiota_arr[j], Z_re[l][m][k][j],
			 Z_re_ipp[l][m][k][j], LOWER_RANGE);
	      Ang_Spline(cosiota_arr[j], Z_im[l][m][k][j],
			 Z_im_ipp[l][m][k][j], LOWER_RANGE);
	    }
	  } // k-loop
	} // m-loop
      } // l-loop
    } else {
      Ang_Spline(cosiota_arr[j], rdot_arr[j], rdot_ipp[j], UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], cosiotadot_arr[j], cosiotadot_ipp[j],
		 UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Edot_arr[j], Edot_ipp[j], UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Lzdot_arr[j], Lzdot_ipp[j], UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Qdot_arr[j], Qdot_ipp[j], UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Omega_th_arr[j], Omega_th_ipp[j],
		 UPPER_RANGE);
      Ang_Spline(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j],
		 UPPER_RANGE);
      for (l = 2; l <= lmax; l++) {
	for (m = -lmax; m <= lmax; m++) {
	  for (k = -kmax; k <= kmax; k++) {
	    Ang_Spline(cosiota_arr[j], Z_re[l][m][k][j],
		       Z_re_ipp[l][m][k][j], UPPER_RANGE);
	    Ang_Spline(cosiota_arr[j], Z_im[l][m][k][j],
		       Z_im_ipp[l][m][k][j], UPPER_RANGE);
	  } // k-loop
	} // m-loop
      } // l-loop
    }
  } // j-loop
}

CID::~CID()
{
  if (lmax > 0) {
    free_Real5tensor(Z_re, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    free_Real5tensor(Z_im, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    free_Real5tensor(Spheroid, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    //
    free_Real5tensor(Z_re_ipp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    free_Real5tensor(Z_im_ipp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    free_Real5tensor(Spheroid_ipp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax,
		     1, 2*imax);
    //
    free_Real4tensor(Z_re_r, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    free_Real4tensor(Z_im_r, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    free_Real4tensor(Spheroid_r, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    //
    free_Real4tensor(Z_re_rpp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    free_Real4tensor(Z_im_rpp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
    free_Real4tensor(Spheroid_rpp, 2, lmax, -lmax, lmax, -kmax, kmax, 1, jmax);
  }
  //
  free_Realvector(Qdot_rpp, 1, jmax);
  free_Realvector(Lzdot_rpp, 1, jmax);
  free_Realvector(Edot_rpp, 1, jmax);
  free_Realvector(cosiotadot_rpp, 1, jmax);
  free_Realvector(rdot_rpp, 1, jmax);
  free_Realvector(Omega_th_rpp, 1, jmax);
  free_Realvector(Omega_ph_rpp, 1, jmax);
  //
  free_Realvector(Qdot_r, 1, jmax);
  free_Realvector(Lzdot_r, 1, jmax);
  free_Realvector(Edot_r, 1, jmax);
  free_Realvector(cosiotadot_r, 1, jmax);
  free_Realvector(rdot_r, 1, jmax);
  free_Realvector(Omega_th_r, 1, jmax);
  free_Realvector(Omega_ph_r, 1, jmax);
  //
  free_Realvector(r_lso_arr, 1, LSO_POINTS);
  free_Realvector(cosiota_lso_arr, 1, LSO_POINTS);
  free_Realvector(r_lso_ipp, 1, LSO_POINTS);
  //
  free_Realmatrix(Qdot_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(Lzdot_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(Edot_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(cosiotadot_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(rdot_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(Omega_th_ipp, 1, jmax, 1, 2*imax);
  free_Realmatrix(Omega_ph_ipp, 1, jmax, 1, 2*imax);
  //
  free_Realmatrix(Qdot_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(Lzdot_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(Edot_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(cosiotadot_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(rdot_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(Omega_th_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(Omega_ph_arr, 1, jmax, 1, 2*imax);
  free_Realmatrix(cosiota_arr, 1, jmax, 1, 2*imax);
  //
  free_Realvector(rstar_arr, 1, jmax);
  free_Realvector(r_arr, 1, jmax);
  //
  free_ivector(ret_exists, 1, jmax);
}

void CID::Spheroids(const Real costheta_view)
{
  if (lmax > 0) {
    for (l = 2; l <= lmax; l++) {
      for (m = -l; m <= l; m++) {
	for (k = -kmax; k <= kmax; k++) {
	  for (j = 1; j <= jmax; j++) {
	    for (i = 1; i <= 2*imax; i++) {
	      SWSH swsh(l, -2, m, a*m*Omega_ph_arr[j][i] +
			a*k*Omega_th_arr[j][i]);
	      Spheroid[l][m][k][j][i] = swsh.spheroid(costheta_view);
	    } // i-loop
	    if (ret_exists[j]) {
	      if (m == 0) {
		Ang_Spline(cosiota_arr[j], Spheroid[l][m][k][j],
			   Spheroid_ipp[l][m][k][j], FULL_RANGE);
	      } else { // m != 0
		//
		// Split into subranges
		//
		Ang_Spline(cosiota_arr[j], Spheroid[l][m][k][j],
			   Spheroid_ipp[l][m][k][j], UPPER_RANGE);
		Ang_Spline(cosiota_arr[j], Spheroid[l][m][k][j],
			   Spheroid_ipp[l][m][k][j], LOWER_RANGE);
	      }
	    } else { // !ret_exists[j]
	      Ang_Spline(cosiota_arr[j], Spheroid[l][m][k][j],
			 Spheroid_ipp[l][m][k][j], UPPER_RANGE);
	    } // else
	  } // j-loop
	} // k-loop
      } // m-loop
    } // l-loop
  } else return;
}

void CID::Pfluxes()
{
//   if (lmax > 0) {
//     for (i = 1; i <= imax; i++) {
//       Complex **Z; Z = Complexmatrix(2, lmax, -lmax, lmax);
//       for (l = 2; l <= lmax; l++)
// 	for (m = -l; m <= l; m++)
// 	  Z[l][m] = Z_re[l][m][i]*exp(II*Z_im[l][m][i]);
//       //
//       Complex tmp_Pplus = Pplus(Z, Omega_arr[i])/(4.*M_PI);
//       Pplus_r[i] = tmp_Pplus.real();
//       Pplus_i[i] = tmp_Pplus.imag();
//       Complex tmp_Pminus = Pminus(Z, Omega_arr[i])/(4.*M_PI);
//       Pminus_r[i] = tmp_Pminus.real();
//       Pminus_i[i] = tmp_Pminus.imag();
//       Complex tmp_Pz = Pz(Z, Omega_arr[i])/(4.*M_PI);
//       Pz_r[i] = tmp_Pz.real();
//       Pz_i[i] = tmp_Pz.imag();
//       //
//       free_Complexmatrix(Z, 2, lmax, -lmax, lmax);
//     }
//     spline(rstar_arr, Pplus_r, imax, 1.e30, 1.e30, Pplus_r_pp);
//     spline(rstar_arr, Pplus_i, imax, 1.e30, 1.e30, Pplus_i_pp);
//     spline(rstar_arr, Pminus_r, imax, 1.e30, 1.e30, Pminus_r_pp);
//     spline(rstar_arr, Pminus_i, imax, 1.e30, 1.e30, Pminus_i_pp);
//     spline(rstar_arr, Pz_r, imax, 1.e30, 1.e30, Pz_r_pp);
//     spline(rstar_arr, Pz_i, imax, 1.e30, 1.e30, Pz_i_pp);
//   } else return;
}

void CID::Get_lso(const Real cosincl)
{
  splint(cosiota_lso_arr, r_lso_arr, r_lso_ipp, LSO_POINTS, cosincl, &r_lso);
}

void CID::Get_coordsdot_omegas(const Real r, const Real cosincl)
{
  const Real rs = Kerr::rstar(r, a);
  rad_offset = 0;
  //
  // At each radius, get data at this inclination.  Only include data
  // inside the LSO radius at this inclination.
  //
  Get_lso(cosincl);
  for (j = 1; j <= jmax; j++) {
    if (r_arr[j] >= r_lso) {
      if (ret_exists[j]) {
	Ang_Splint(cosiota_arr[j], rdot_arr[j], rdot_ipp[j], cosincl,
		   rdot_r[j], FULL_RANGE);
	Ang_Splint(cosiota_arr[j], cosiotadot_arr[j], cosiotadot_ipp[j],
		   cosincl, cosiotadot_r[j], FULL_RANGE);
	Ang_Splint(cosiota_arr[j], Omega_th_arr[j], Omega_th_ipp[j], cosincl,
		   Omega_th_r[j], FULL_RANGE);
	//
	// Omega_ph has subranges
	//
	if (cosincl > 0)
	  Ang_Splint(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j], cosincl,
		     Omega_ph_r[j], UPPER_RANGE);
	else
	  Ang_Splint(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j], cosincl,
		     Omega_ph_r[j], LOWER_RANGE);
      } else { // !ret_exists[j]
	Ang_Splint(cosiota_arr[j], rdot_arr[j], rdot_ipp[j], cosincl,
		   rdot_r[j], UPPER_RANGE);
	Ang_Splint(cosiota_arr[j], cosiotadot_arr[j], cosiotadot_ipp[j],
		   cosincl, cosiotadot_r[j], UPPER_RANGE);
	Ang_Splint(cosiota_arr[j], Omega_th_arr[j], Omega_th_ipp[j], cosincl,
		   Omega_th_r[j], UPPER_RANGE);
	Ang_Splint(cosiota_arr[j], Omega_ph_arr[j], Omega_ph_ipp[j], cosincl,
		   Omega_ph_r[j], UPPER_RANGE);
      }
      if (fabs(cosiotadot_r[j]) < 1.e-14) cosiotadot_r[j] = 0.;
    } else rad_offset++;
  }
  //
  // Compute fields at current point.
  //
  Rad_Spline(rdot_r, rdot_rpp, rad_offset);
  Rad_Splint(rdot_r, rdot_rpp, rs, rdot, rad_offset);
  rdot /= (r - r_lso);
  //
  Rad_Spline(cosiotadot_r, cosiotadot_rpp, rad_offset);
  Rad_Splint(cosiotadot_r, cosiotadot_rpp, rs, cosiotadot, rad_offset);
  //
  Rad_Spline(Omega_ph_r, Omega_ph_rpp, rad_offset);
  Rad_Splint(Omega_ph_r, Omega_ph_rpp, rs, Omega_ph, rad_offset);
  //
  Rad_Spline(Omega_th_r, Omega_th_rpp, rad_offset);
  Rad_Splint(Omega_th_r, Omega_th_rpp, rs, Omega_th, rad_offset);
}

void CID::Get_fluxes(const Real r, const Real cosincl)
{
  const Real rs = Kerr::rstar(r, a);
  rad_offset = 0;
  //
  // At each radius, get data at this inclination.  Only include data
  // inside the LSO radius at this inclination.
  //
  Get_lso(cosincl);
  for (j = 1; j <= jmax; j++) {
    if (r_arr[j] >= r_lso) {
      if (ret_exists[j]) {
	Ang_Splint(cosiota_arr[j], Edot_arr[j], Edot_ipp[j], cosincl,
		   Edot_r[j], FULL_RANGE);
	Ang_Splint(cosiota_arr[j], Lzdot_arr[j], Lzdot_ipp[j], cosincl,
		   Lzdot_r[j], FULL_RANGE);
	Ang_Splint(cosiota_arr[j], Qdot_arr[j], Qdot_ipp[j], cosincl,
		   Qdot_r[j], FULL_RANGE);
      } else {
	Ang_Splint(cosiota_arr[j], Edot_arr[j], Edot_ipp[j], cosincl,
		   Edot_r[j], UPPER_RANGE);
	Ang_Splint(cosiota_arr[j], Lzdot_arr[j], Lzdot_ipp[j], cosincl,
		   Lzdot_r[j], UPPER_RANGE);
	Ang_Splint(cosiota_arr[j], Qdot_arr[j], Qdot_ipp[j], cosincl,
		   Qdot_r[j], UPPER_RANGE);
      }
    } else rad_offset++;
  }
  //
  // Compute fields at current point.
  //
  Rad_Spline(Edot_r, Edot_rpp, rad_offset);
  Rad_Splint(Edot_r, Edot_rpp, rs, Edot, rad_offset);
  Rad_Spline(Lzdot_r, Lzdot_rpp, rad_offset);
  Rad_Splint(Lzdot_r, Lzdot_rpp, rs, Lzdot, rad_offset);
  Rad_Spline(Qdot_r, Qdot_rpp, rad_offset);
  Rad_Splint(Qdot_r, Qdot_rpp, rs, Qdot, rad_offset);
}

void CID::Get_wave(const Real r, const Real cosincl,
		   const Real Nph, const Real Nth, const Real phi_view)
{
  const Real rs = Kerr::rstar(r, a);
  //
  // At each radius, get data at this inclination.  Only include data
  // inside the LSO radius at this inclination.
  //
  Get_lso(cosincl);
  hp = 0.; hc = 0.;
  Real Z_re_int, Z_im_int, Spheroid_int, dhc, dhp;
  Complex Zhere;
  RRGW rrgw;
  //
  for (l = 2; l <= lmax; l++) {
    for (m = -l; m <= l; m++) {
      for (k = -kmax; k <= kmax; k++) {
	if (m != 0 || k != 0) {
	  rad_offset = 0;
	  for (j = 1; j <= jmax; j++) {
	    if (r_arr[j] >= r_lso) {
	      if (ret_exists[j]) {
		if (m == 0) {
		  Ang_Splint(cosiota_arr[j], Z_re[l][m][k][j],
			     Z_re_ipp[l][m][k][j], cosincl,
			     Z_re_r[l][m][k][j], FULL_RANGE);
		  Ang_Splint(cosiota_arr[j], Z_im[l][m][k][j],
			     Z_im_ipp[l][m][k][j], cosincl,
			     Z_im_r[l][m][k][j], FULL_RANGE);
		  Ang_Splint(cosiota_arr[j], Spheroid[l][m][k][j],
			     Spheroid_ipp[l][m][k][j], cosincl,
			     Spheroid_r[l][m][k][j], FULL_RANGE);
 		} else { // m != 0
		  //
		  // Need to handle subranges
		  //
		  if (cosincl > 0) {
		    Ang_Splint(cosiota_arr[j], Z_re[l][m][k][j],
			       Z_re_ipp[l][m][k][j], cosincl,
			       Z_re_r[l][m][k][j], UPPER_RANGE);
		    Ang_Splint(cosiota_arr[j], Z_im[l][m][k][j],
			       Z_im_ipp[l][m][k][j], cosincl,
			       Z_im_r[l][m][k][j], UPPER_RANGE);
		    Ang_Splint(cosiota_arr[j], Spheroid[l][m][k][j],
			       Spheroid_ipp[l][m][k][j], cosincl,
			       Spheroid_r[l][m][k][j], UPPER_RANGE);
		  } else { // cosincl < 0
		    Ang_Splint(cosiota_arr[j], Z_re[l][m][k][j],
			       Z_re_ipp[l][m][k][j], cosincl,
			       Z_re_r[l][m][k][j], LOWER_RANGE);
		    Ang_Splint(cosiota_arr[j], Z_im[l][m][k][j],
			       Z_im_ipp[l][m][k][j], cosincl,
			       Z_im_r[l][m][k][j], LOWER_RANGE);
		    Ang_Splint(cosiota_arr[j], Spheroid[l][m][k][j],
			       Spheroid_ipp[l][m][k][j], cosincl,
			       Spheroid_r[l][m][k][j], LOWER_RANGE);
		  }
		}
	      } else { // !ret_exists[j]
		Ang_Splint(cosiota_arr[j], Z_re[l][m][k][j],
			   Z_re_ipp[l][m][k][j], cosincl,
			   Z_re_r[l][m][k][j], UPPER_RANGE);
		Ang_Splint(cosiota_arr[j], Z_im[l][m][k][j],
			   Z_im_ipp[l][m][k][j], cosincl,
			   Z_im_r[l][m][k][j], UPPER_RANGE);
		Ang_Splint(cosiota_arr[j], Spheroid[l][m][k][j],
			   Spheroid_ipp[l][m][k][j], cosincl,
			   Spheroid_r[l][m][k][j], UPPER_RANGE);
	      }
	    } else rad_offset++;
	  } // j-loop
	  //
	  // Compute fields at current point.
	  //
	  Rad_Spline(Z_re_r[l][m][k], Z_re_rpp[l][m][k], rad_offset);
	  Rad_Splint(Z_re_r[l][m][k], Z_re_rpp[l][m][k], rs, Z_re_int,
		     rad_offset);
	  Rad_Spline(Z_im_r[l][m][k], Z_im_rpp[l][m][k], rad_offset);
	  Rad_Splint(Z_im_r[l][m][k], Z_im_rpp[l][m][k], rs, Z_im_int,
		     rad_offset);
	  Rad_Spline(Spheroid_r[l][m][k], Spheroid_rpp[l][m][k], rad_offset);
	  Rad_Splint(Spheroid_r[l][m][k], Spheroid_rpp[l][m][k], rs,
		     Spheroid_int, rad_offset);
	  Zhere = Z_re_int + II*Z_im_int;
	  rrgw.Wave(m, k, Nph, Nth, phi_view, Spheroid_int,
		    m*Omega_ph + k*Omega_th, Zhere, dhp, dhc);
	  hp += dhp; hc += dhc;
	} // if (m != 0 && k != 0)
      } // k-loop
    } // m-loop
  } // l-loop
}

void CID::Ang_Spline(Real *x, Real *y, Real *y_xpp, const int Type)
{
  Real deriv_lo, deriv_hi;
  if (Type == FULL_RANGE) {
    //
    // If low end touches equator, set deriv_lo to 0.  Otherwise, set
    // to 1.e30 (natural spline).  High end always touches equator.
    //
    if (x[1] == -1.) deriv_lo = 0.;
    else deriv_lo = 1.e30;
    deriv_hi = 0.;
    //
    spline(x, y, 2*imax, deriv_lo, deriv_hi, y_xpp);
  } else if (Type == LOWER_RANGE) {
    //
    // If low end touches equator, set deriv_lo to 0.  Otherwise, set
    // to 1.e30 (natural spline).  High end never touches equator.
    //
    if (x[1] == -1.) deriv_lo = 0.;
    else deriv_lo = 1.e30;
    deriv_hi = 1.e30;
    //
    spline(x, y, imax, deriv_lo, deriv_hi, y_xpp);
  } else if (Type == UPPER_RANGE) {
    //
    // Low end never touches equator.  High end always touches
    // equator.
    //
    deriv_lo = 1.e30;
    deriv_hi = 0.;
    spline(x + imax, y + imax, imax, deriv_lo, deriv_hi, y_xpp + imax);
  } else
    Die("Unexpected range type in Ang_Spline!!");
}

void CID::Rad_Spline(Real *y, Real *y_rpp, const int rad_offset)
{
  spline(rstar_arr + rad_offset, y + rad_offset, jmax - rad_offset,
	 1.e30, 1.e30, y_rpp + rad_offset);
}

void CID::Ang_Splint(Real *x, Real *y, Real *y_xpp, const Real xc,
		     Real & yval, const int Type)
{
  if (Type == FULL_RANGE) {
    splint(x, y, y_xpp, 2*imax, xc, &yval);
  } else if (Type == LOWER_RANGE) {
    splint(x, y, y_xpp, imax, xc, &yval);
  } else if (Type == UPPER_RANGE) {
    splint(x + imax, y + imax, y_xpp + imax, imax, xc, &yval);
  } else
    Die("Unexpected range type in Ang_Splint!!");
}

void CID::Rad_Splint(Real *y, Real *y_rpp, const Real rs, Real & yval,
		     const int rad_offset)
{
  splint(rstar_arr + rad_offset, y + rad_offset, y_rpp + rad_offset,
	 jmax - rad_offset, rs, &yval);
}

// // This code will not get your mother.  Sorry.
// void CID::Get_mom(const Real r, const Real phi)
// {
//   const Real rs = Kerr::rstar(r, a);

//   Real Pplus_r_int, Pplus_i_int;
//   Real Pminus_r_int, Pminus_i_int;
//   Real Pz_r_int, Pz_i_int;
//   //
//   splint(rstar_arr, Pplus_r, Pplus_r_pp, imax, rs, &Pplus_r_int);
//   splint(rstar_arr, Pplus_i, Pplus_i_pp, imax, rs, &Pplus_i_int);
//   splint(rstar_arr, Pminus_r, Pminus_r_pp, imax, rs, &Pminus_r_int);
//   splint(rstar_arr, Pminus_i, Pminus_i_pp, imax, rs, &Pminus_i_int);
//   splint(rstar_arr, Pz_r, Pz_r_pp, imax, rs, &Pz_r_int);
//   splint(rstar_arr, Pz_i, Pz_i_pp, imax, rs, &Pz_i_int);
//   //
//   Complex plus = exp(II*phi)*(Pplus_r_int + II*Pplus_i_int);
//   Complex minus = exp(-II*phi)*(Pminus_r_int + II*Pminus_i_int);
//   //
//   // The next three combinations should all actually be _real_
//   // numbers.  We leave them like this for now to make sure nothing is
//   // acting screwy.  Note the sign --- we want system's recoil,
//   // opposite of the momentum carried by the radiation.
//   //
//   pdot_x = -0.5*(plus + minus);
//   pdot_y = II*0.5*(plus - minus);
//   pdot_z = -(Pz_r_int + II*Pz_i_int);
// }

// Complex CID::Pplus(Complex **Z, const Real omega)
// {
//   int l, lp, m, mp;
//   Complex sum = 0.;
//   for (l = 2; l <= lmax; l++) {
//     for (lp = 2; lp <= lmax; lp++) {
//       for (m = -l; m <= l; m++) {
// 	mp = m + 1;
// 	//
// 	// Note the values we exclude.  This is because
// 	// ZH_{lm} = 0 for m = 0
// 	// |m + 1| <= lp is required
// 	//
//  	if (m != 0 && mp != 0 && abs(mp) <= lp) {
// 	  Complex Zed_l = Z[l][m];
// 	  Complex Zed_lp = Z[lp][m + 1];
	  
// 	  sum += Zed_l*conj(Zed_lp)*Il_lp_m_mp(l, lp, m, m + 1, omega)/
// 	    (m*(m + 1)*omega*omega);
// 	}
//       }
//     }
//   }
//   return(sum);
// }

// Complex CID::Pminus(Complex **Z, const Real omega)
// {
//   int l, lp, m, mp;
//   Complex sum = 0.;
//   for (l = 2; l <= lmax; l++) {
//     for (lp = 2; lp <= lmax; lp++) {
//       for (m = -l; m <= l; m++) {
// 	mp = m - 1;
// 	//
// 	// Note the values we exclude.  This is because
// 	// ZH_{lm} = 0 for m = 0
// 	// |m - 1| <= lp is required
// 	//
//  	if (m != 0 && mp != 0 && abs(mp) <= lp) {
// 	  Complex Zed_l = Z[l][m];
// 	  Complex Zed_lp = Z[lp][m - 1];

// 	  sum += Zed_l*conj(Zed_lp)*Il_lp_m_mp(l, lp, m, m - 1, omega)/
// 	    (m*(m - 1)*omega*omega);
//  	}
//       }
//     }
//   }
//   return(sum);
// }

// Complex CID::Pz(Complex **Z, const Real omega)
// {
//   int l, lp, m;
//   Complex sum = 0;
//   for (l = 2; l <= lmax; l++) {
//     for (lp = 2; lp <= lmax; lp++) {
//       for (m = -l; m <= l; m++) {
// 	//
// 	// Note the values we exclude.  This is because
// 	// ZH_{lm} = 0 for m = 0
// 	// |m| <= lp is required
// 	//
//  	if (m != 0 && abs(m) <= lp) {
// 	  Complex Zed_l = Z[l][m];
// 	  Complex Zed_lp = Z[lp][m];

// 	  sum += Zed_l*conj(Zed_lp)*Jl_lp_m(l, lp, m, omega)/
// 	    (m*m*omega*omega);
//  	}
//       }
//     }
//   }
//   if (abs(sum) < 1.e-15) sum = 0.;
//   return(sum);
// }

// Real CID::Il_lp_m_mp(const int l, const int lp, const int m, const int mp,
// 		      const Real omega)
// {
//   const Real aw = m*a*omega;
//   const Real awp = mp*a*omega;
//   //
//   SWSH lm(l, -2, m, aw);
//   SWSH lpmp(lp, -2, mp, awp);
//   Clebsch cg;
//   Real sum = 0.;
//   int i, ip;
//   int k, kp;
//   //
//   // Note some slightly complicated shenanigans required because the
//   // way that the coefficients are stored.  b[1] corresponds to the
//   // entry for lmin, b[2] corresponds to lmin + 1, etc.
//   //
//   for (i = 1; i <= lm.N; i++) {
//     k = i - 1 + lm.lmin;
//     //
//     // First term
//     //
//     kp = k - 1;
//     // Find corresponding index.
//     ip = kp + 1 - lpmp.lmin;
//     if (ip > 0 && ip <= lpmp.N) {
//       sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
//     }
//     //
//     // Second term
//     //
//     kp = k;
//     ip = kp + 1 - lpmp.lmin;
//     if (ip > 0 && ip <= lpmp.N) {
//       sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
//     }
//     //
//     // Third term
//     //
//     kp = k + 1;
//     ip = kp + 1 - lpmp.lmin;
//     if (ip > 0 && ip <= lpmp.N) {
//       sum += lm.b[i]*lpmp.b[ip]*cg.sinthetabrac(-2, k, kp, m, mp);
//     }
//   }
//   return(sum);
// }

// Real CID::Jl_lp_m(const int l, const int lp, const int m, const Real omega)
// {
//   const Real aw = m*a*omega;
//   //
//   SWSH lm(l, -2, m, aw);
//   SWSH lpm(lp, -2, m, aw);
//   Clebsch cg;
//   Real sum = 0.;
//   int i, ip;
//   int k, kp;
//   //
//   // Note some slightly complicated shenanigans required because the
//   // way that the coefficients are stored.  b[1] corresponds to the
//   // entry for lmin, b[2] corresponds to lmin + 1, etc.
//   //
//   for (i = 1; i <= lm.N; i++) {
//     k = i - 1 + lm.lmin;
//     //
//     // First term
//     //
//     kp = k - 1;
//     // Find corresponding index.
//     ip = kp + 1 - lpm.lmin;
//     if (ip > 0 && ip <= lpm.N) {
//       sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
//     }
//     //
//     // Second term
//     //
//     kp = k;
//     ip = kp + 1 - lpm.lmin;
//     if (ip > 0 && ip <= lpm.N) {
//       sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
//     }
//     //
//     // Third term
//     //
//     kp = k + 1;
//     ip = kp + 1 - lpm.lmin;
//     if (ip > 0 && ip <= lpm.N) {
//       sum += lm.b[i]*lpm.b[ip]*cg.xbrac(-2, k, kp, m);
//     }
//   }
//   return(sum);
// }
