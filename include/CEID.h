//---------------------------------------------------------------------------
//
// $Id: CEID.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Circular, equatorial inspiral data.
//
#ifndef _CEID_H
#define _CEID_H

#include <math.h>
#include "Globals.h"
#include "SWSH.h"
#include "RRGW.h"

void spline(Real xa[], Real fa[], int n, Real fp1, Real fpn, Real fpp[]);
void splint(Real xa[], Real fa[], Real fpp[], int n, Real x, Real *f);

class CEID {
public:
  CEID(char *basename, const Real arrmin, const Real rmax, const Real dr,
       const int ellmax);
  ~CEID();

  Real rdot, Omega;
  Real Omega_max, rmin;
  Real Edot, Lzdot;
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
  // Get stuff, interpolated to current coordinate in parameter space
  //
  void Get_rdot_omega(const Real r);
  void Get_fluxes(const Real r);
  void Get_wave(const Real r, const Real phi, const Real phi_view);
  void Get_mom(const Real r, const Real phi);

private:
  void Allocate();
  void ReadIn(char *inname);
  void Make_radial_splines();

  Real a, r_isco;
  int lmax;
  int j, jmax, l, m;
  //
  // {stuff}_pp[] are the second derivative arrays created by spline
  // and then used by splint.
  //
  Real *tmp_r, *tmp_Omega, *tmp_rdot, *tmp_Edot, *tmp_Lzdot;
  Real *r_arr, *rstar_arr, *Omega_arr, *rdot_arr, *Edot_arr, *Lzdot_arr;
  Real *Omega_pp, *rdot_pp, *Edot_pp, *Lzdot_pp;
  //
  // Used for waveforms.
  //
  Real ***tmp_Z_re, ***tmp_Z_im;
  Real ***Z_re, ***Z_im, ***Spheroid;
  Real ***Z_re_pp, ***Z_im_pp, ***Spheroid_pp;
  //
  // Used for recoil.
  //
  Real *Pplus_r, *Pminus_r, *Pz_r;
  Real *Pplus_i, *Pminus_i, *Pz_i;
  Real *Pplus_r_pp, *Pminus_r_pp, *Pz_r_pp;
  Real *Pplus_i_pp, *Pminus_i_pp, *Pz_i_pp;
  //
  // Functions needed to compute momentum flux.
  //
  Complex Pplus(Complex **Z, const Real omega);
  Complex Pminus(Complex **Z, const Real omega);
  Complex Pz(Complex **Z, const Real omega);
  Real Il_lp_m_mp(const int l, const int lp, const int m, const int mp,
		  const Real omega);
  Real Jl_lp_m(const int l, const int lp, const int m, const Real omega);
};
#endif
