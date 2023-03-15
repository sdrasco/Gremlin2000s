//---------------------------------------------------------------------------
//
// $Id: CID.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial inspiral data.
//
#ifndef _CID_H
#define _CID_H

#include <math.h>
#include "Globals.h"
#include "CKG.h"
#include "LB.h"
#include "SWSH.h"
#include "RRGW.h"

//
// Number of points used to define r_LSO(cosiota)
//
#define LSO_POINTS 1000
//
// Useful defines/flags
//
#define PROGRADE 1
#define RETROGRADE 0
//
// Three spline types:
// FULL_RANGE: data ranges from 1 ... 2 imax
// LOWER_RANGE: data ranges from 1 ... imax
// UPPER_RANGE: data ranges from imax + 1 ... 2 imax
//
#define FULL_RANGE 1
#define LOWER_RANGE 2
#define UPPER_RANGE 3

void spline(Real xa[], Real fa[], int n, Real fp1, Real fpn, Real fpp[]);
void splint(Real xa[], Real fa[], Real fpp[], int n, Real x, Real *f);

class CID {
public:
  CID(char *basename, const int ayemax, const Real arrmin, const Real rmax,
      const Real dr, const int ellmax, const int kaymax);
  ~CID();
  //
  // The following two things must be public in order to check whether
  // there is enough data to do a run.
  //
  int jmax;
  Real *r_arr;
  //
  Real a;
  Real Omega_max, rmin;
  Real r_lso;
  //
  Real rdot, cosiotadot, Omega_ph, Omega_th;
  Real Edot, Lzdot, Qdot;
  Real hp, hc;
  Complex pdot_x, pdot_y, pdot_z;
  //
  // Loads spheroid data.
  //
  void Spheroids(const Real costheta_view);
  //
  // Loads momentum flux data.
  //
  void Pfluxes();
  //
  // Get the inclination of the LSO at this radius
  //
  void Get_lso(const Real incl);
  //
  // Get stuff, interpolated to current coordinate in parameter space
  //
  void Get_coordsdot_omegas(const Real r, const Real cosincl);
  void Get_fluxes(const Real r, const Real cosincl);
  void Get_wave(const Real r, const Real cosincl,
		const Real Nph, const Real Nth, const Real phi_view);
  //   void Get_mom(const Real r, const Real phi);

private:
  void Allocate();
  void ReadIn(char *inname, const int proflag);
  void Make_lso();
  void Make_angular_splines();
  //
  // Wrapper for splining angular data
  //
  void Ang_Spline(Real *x, Real *y, Real *y_xpp, const int Type);
  void Rad_Spline(Real *y, Real *y_rpp, const int rad_offset);
  void Ang_Splint(Real *x, Real *y, Real *y_xpp, const Real xc, Real & yval,
		  const int Type);
  void Rad_Splint(Real *y, Real *y_rpp, const Real rs, Real & yval,
		  const int rad_offset);

  int i, imax, j, l, m, k;
  int rad_offset;
  int lmax, kmax;
  //
  int *ret_exists;
  //
  Real *rstar_arr;
  Real **cosiota_arr;
  //
  // {stuff}_ipp[] are the second derivative in cosiota arrays created
  // by spline and then used by splint.
  //  
  Real **Omega_ph_arr, **Omega_th_arr;
  Real **rdot_arr, **cosiotadot_arr;
  Real **Edot_arr, **Lzdot_arr, **Qdot_arr;
  Real **Omega_ph_ipp, **Omega_th_ipp;
  Real **rdot_ipp, **cosiotadot_ipp;
  Real **Edot_ipp, **Lzdot_ipp, **Qdot_ipp;
  //
  // Radius of LSO depends only on cosiota (index i).
  //
  Real *r_lso_arr, *cosiota_lso_arr;
  Real *r_lso_ipp;
  //
  // All of these data interpolated onto some cosiota, and 2nd derivs
  // in r.
  //
  Real *cosiota_r, *Omega_ph_r, *Omega_th_r, *rdot_r, *cosiotadot_r;
  Real *Edot_r, *Lzdot_r, *Qdot_r;
  Real *Omega_ph_rpp, *Omega_th_rpp, *rdot_rpp, *cosiotadot_rpp;
  Real *Edot_rpp, *Lzdot_rpp, *Qdot_rpp;
  //
  // Used for waveforms.
  //
  Real *****Z_re, *****Z_im, *****Spheroid;
  Real *****Z_re_ipp, *****Z_im_ipp, *****Spheroid_ipp;
  Real ****Z_re_r, ****Z_im_r, ****Spheroid_r;
  Real ****Z_re_rpp, ****Z_im_rpp, ****Spheroid_rpp;
/*   // */
/*   // Used for recoil. */
/*   // */
/*   Real *Pplus_r, *Pminus_r, *Pz_r; */
/*   Real *Pplus_i, *Pminus_i, *Pz_i; */
/*   Real *Pplus_r_pp, *Pminus_r_pp, *Pz_r_pp; */
/*   Real *Pplus_i_pp, *Pminus_i_pp, *Pz_i_pp; */
/*   // */
/*   // Functions needed to compute momentum flux. */
/*   // */
/*   Complex Pplus(Complex **Z, const Real omega); */
/*   Complex Pminus(Complex **Z, const Real omega); */
/*   Complex Pz(Complex **Z, const Real omega); */
/*   Real Il_lp_m_mp(const int l, const int lp, const int m, const int mp, */
/* 		  const Real omega); */
/*   Real Jl_lp_m(const int l, const int lp, const int m, const Real omega); */
};
#endif
