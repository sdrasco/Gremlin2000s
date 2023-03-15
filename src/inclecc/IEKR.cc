//
// Methods of class IEKR (Inclined Eccentric Kerr Radiation)
//
// Computes the expansion coefficients Z_{lmkn} for the inhomogenous 
// part of the radial Teukolsky function.  The source term is an 
// arbitrary bound geodesic (iekg object).
//
// Steve Drasco
//
#include <math.h>
#include "Globals.h"
#include "IEKR.h"
#include "NRUtil.h"
#include "SWSH.h"
#include "IEKG.h"
#include "SNT.h"

// Constructor.
IEKR::IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index, 
          int n_index, const Real tol) : 
  iekg(iekg_in), l(l_index), m(m_index), k(k_index), n(n_index), eps(tol)
{
  // TEST
  cerr << "this is output from IEKR.cc\n";

  // make sure indices are ok
  if(l < 2) Die("Illegal value of l in IEKR");
  if(abs(m) > l) Die("Illegal value of m in IEKR");

  // local storage of stuff from the geodesic object
  a = iekg->a;
  E = iekg->E;
  Lz = iekg->Lz;
  Q = iekg->Q;
  omega = iekg->Omega(m,k,n);
  Gamma = iekg->Gamma;

  // in this constructor, we don't use MaxEI and MaxEH
  MaxEI = MaxEH = 0.0;

  // construct the (-2) spin weighted spheroidal harmonic
  SWSH swsh_tmp(l, -2, m, a*omega); 
  swsh = &swsh_tmp;

  // initialize things for the homogeneous solutions
  // 
  // The number Nspline is the number of points at which the homogenous
  // fields will be computed explicitly using the integrator.  Later, when
  // doing the psi-integral, these are used for a spline interpolation to 
  // give the fields at any point.  
  //
  // Nspline = 25000 is HUGE.  Spline prep takes about 30 seconds for this.
  // Nspline = 5000 gives Z_{lmkn} that agree with 25000 to about 5-9 digits.
  // 
  // In some cases we can get away with far smaller numbers Nspline ~ 10. 
  // 
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps, 1e-10));
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, eps/100.0);
  SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps/100.0,1e-7));
  snt = &snt_init;
  //int Nspline = 30000; 
  //snt->PrepSpline(Nspline);

  // get Almkn is as in (4.10) of PRD 61 084004 (including its sign error)
  Almkn = - 2.0 * II * omega * snt->Bin;

  // set integration boundaries
  IntToPi = 1;

  // Turn integrand scaling on or off.  If on, the radial integrand is rescaled
  // by a change of varriables according to: 
  //
  //  int_0^pi dpsi f(psi)  -->  int_A^B dy exp(-psi/ScaleLambda) f(psi)
  // 
  if(IntToPi)
    ScaleIntegrand = 1;
  else
    ScaleIntegrand = 0;

  // compute alpha_{lmkn} addapted from RRGW::alpha_func([stuff])
  //
  // Factor that appears in the fluxes on the horizon; arises from
  // the conversion of the Newman-Penrose tetrad to a Hartle-
  // Hawking tetrad (used to describe how the horizon grows due to
  // stuff falling in).  See Teukolsky and Press 1974.
  //
  const Real lamb = swsh->lambda;
  const Real w = omega;
  const Real r_plus = Kerr::rplus(a);
  const Real P = w - Real(m)*a/(2.0*r_plus);
  const Real aeps = sqrt(1. - a*a)/(4.*r_plus);
  const Real lp2 = lamb + 2.;
  const Real aw_m_m = a*w - ((Real)m);
  const Real tmp1 = lp2*lp2 - 4.*a*w*aw_m_m;
  const Real tmp2 = lamb*lamb - 36.*a*w*aw_m_m;
  const Real tmp3 = 48.*a*w*(2.*lamb + 3.)*(2.*a*w - ((Real)m));
  const Real tmp4 = 256.*pow(2.*r_plus, 5.)*P*(P*P + 4.*aeps*aeps)*
    (P*P + 16.*aeps*aeps)*w*w*w;
  const Real abssqrClm = tmp1*tmp2 + tmp3 + 144.*w*w*(1. - a*a);
  alpha = tmp4/abssqrClm;

  // determine rescaling rules for the Z^H_{lmkn} integrand
  DoH = 1;
  EstimateScaleLambda(); 
  Real A = 0.0;
  Real B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

  // compute Z^H_{lmkn}
  int PsiIntegral=1;
  if(ScaleIntegrand) {
    ZH = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, eps); 
  } else {
    if(IntToPi)
      ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, eps);
    else
      ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, eps);
  }
  ZH /= -Almkn;
  ZH *= 1.0 / (2.0 * M_PI * Gamma);

  // determine rescaling rules for the Z^inf_{lmkn} integrand
  DoH = 0;
  EstimateScaleLambda();
  B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

  // compute Z^{\infty}_{lmkn}
  if(ScaleIntegrand) {
    ZInf = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, eps);
  } else {
    if(IntToPi) 
      ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, eps); 
    else
      ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, eps);
  }
  ZInf /= -Almkn;
  ZInf *= 1.0 / (2.0 * M_PI * Gamma);

}
// overloaded constructor for passing in MaxEI and MaxEH
IEKR::IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index,
          int n_index, const Real tol, Real MaxEIin, Real MaxEHin) :
  iekg(iekg_in), l(l_index), m(m_index), k(k_index), n(n_index), eps(tol), MaxEI(MaxEIin), MaxEH(MaxEHin)
{
  // make sure indices are ok
  if(l < 2) Die("Illegal value of l in IEKR");
  if(abs(m) > l) Die("Illegal value of m in IEKR");

  // local storage of stuff from the geodesic object
  a = iekg->a;
  E = iekg->E;
  Lz = iekg->Lz;
  Q = iekg->Q;
  omega = iekg->Omega(m,k,n);
  Gamma = iekg->Gamma;

  // construct the (-2) spin weighted spheroidal harmonic
  SWSH swsh_tmp(l, -2, m, a*omega);
  swsh = &swsh_tmp;

  // initialize things for the homogeneous solutions
  //
  // The number Nspline is the number of points at which the homogenous
  // fields will be computed explicitly using the integrator.  Later, when
  // doing the psi-integral, these are used for a spline interpolation to
  // give the fields at any point.
  //
  // Nspline = 25000 is HUGE.  Spline prep takes about 30 seconds for this.
  // Nspline = 5000 gives Z_{lmkn} that agree with 25000 to about 5-9 digits.
  //
  // In some cases we can get away with far smaller numbers Nspline ~ 10.
  //
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps, 1e-10));
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, eps/100.0);
  SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps/100.0,1e-7));
  snt = &snt_init;
  //int Nspline = 30000;
  //snt->PrepSpline(Nspline);

  // get Almkn is as in (4.10) of PRD 61 084004 (including its sign error)
  Almkn = - 2.0 * II * omega * snt->Bin;

  // set integration boundaries
  IntToPi = 1;

  // Turn integrand scaling on or off.  If on, the radial integrand is rescaled
  // by a change of varriables according to:
  //
  //  int_0^pi dpsi f(psi)  -->  int_A^B dy exp(-psi/ScaleLambda) f(psi)
  //
  if(IntToPi)
    ScaleIntegrand = 1;
  else
    ScaleIntegrand = 0;

  // compute alpha_{lmkn} addapted from RRGW::alpha_func([stuff])
  //
  // Factor that appears in the fluxes on the horizon; arises from
  // the conversion of the Newman-Penrose tetrad to a Hartle-
  // Hawking tetrad (used to describe how the horizon grows due to
  // stuff falling in).  See Teukolsky and Press 1974.
  //
  const Real lamb = swsh->lambda;
  const Real w = omega;
  const Real r_plus = Kerr::rplus(a);
  const Real P = w - Real(m)*a/(2.0*r_plus);
  const Real aeps = sqrt(1. - a*a)/(4.*r_plus);
  const Real lp2 = lamb + 2.;
  const Real aw_m_m = a*w - ((Real)m);
  const Real tmp1 = lp2*lp2 - 4.*a*w*aw_m_m;
  const Real tmp2 = lamb*lamb - 36.*a*w*aw_m_m;
  const Real tmp3 = 48.*a*w*(2.*lamb + 3.)*(2.*a*w - ((Real)m));
  const Real tmp4 = 256.*pow(2.*r_plus, 5.)*P*(P*P + 4.*aeps*aeps)*
    (P*P + 16.*aeps*aeps)*w*w*w;
  const Real abssqrClm = tmp1*tmp2 + tmp3 + 144.*w*w*(1. - a*a);
  alpha = tmp4/abssqrClm;

  // determine rescaling rules for the Z^H_{lmkn} integrand
  DoH = 1;
  EstimateScaleLambda();
  Real A = 0.0;
  Real B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

  // compute Z^H_{lmkn}
  int PsiIntegral=1;
  if(ScaleIntegrand) {
    ZH = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, Min(1e-5,eps));
  } else {
    if(IntToPi)
      ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, Min(1e-5,eps));
    else
      ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, Min(1e-5,eps));
  }
  ZH /= -Almkn;
  ZH *= 1.0 / (2.0 * M_PI * Gamma);

  // determine rescaling rules for the Z^inf_{lmkn} integrand
  DoH = 0;
  EstimateScaleLambda();
  B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

  // compute Z^{\infty}_{lmkn}
  if(ScaleIntegrand) {
    ZInf = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, Min(1e-5,eps));
  } else {
    if(IntToPi)
      ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, Min(1e-5,eps));
    else
      ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, Min(1e-5,eps));
  }
  ZInf /= -Almkn;
  ZInf *= 1.0 / (2.0 * M_PI * Gamma);

}

// overloaded constructor: allows calculation of only Z^H_{lmkn} or only Z^{\infty}_{lmkn}
//                         allows for tracking MaxEH and MaxEI
IEKR::IEKR(IEKG *iekg_in, int l_index, int m_index, int k_index,
          int n_index, const Real tol, Real MaxEIin, Real MaxEHin, int DoZH, int DoZI) :
  iekg(iekg_in), l(l_index), m(m_index), k(k_index), n(n_index), eps(tol), MaxEI(MaxEIin), MaxEH(MaxEHin)
{
  // make sure indices are ok
  if(l < 2) Die("Illegal value of l in IEKR");
  if(abs(m) > l) Die("Illegal value of m in IEKR");

  // local storage of stuff from the geodesic object
  a = iekg->a;
  E = iekg->E;
  Lz = iekg->Lz;
  Q = iekg->Q;
  omega = iekg->Omega(m,k,n);
  Gamma = iekg->Gamma;

  // construct the (-2) spin weighted spheroidal harmonic
  SWSH swsh_tmp(l, -2, m, a*omega);
  swsh = &swsh_tmp;

  // initialize things for the homogeneous solutions
  //
  // The number Nspline is the number of points at which the homogenous
  // fields will be computed explicitly using the integrator.  Later, when
  // doing the psi-integral, these are used for a spline interpolation to
  // give the fields at any point.
  //
  // Nspline = 25000 is HUGE.  Spline prep takes about 30 seconds for this.
  // Nspline = 5000 gives Z_{lmkn} that agree with 25000 to about 5-9 digits.
  //
  // In some cases we can get away with far smaller numbers Nspline ~ 10.
  //
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps, 1e-10));
  //SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, eps/100.0);
  SNT snt_init(m, iekg->p, iekg->e, a, omega, swsh->lambda, Min(eps/100.0,1e-7));
  snt = &snt_init;
  //int Nspline = 30000;
  //snt->PrepSpline(Nspline);

  // record Ain
  Ain = snt->Ain;

  // get Almkn is as in (4.10) of PRD 61 084004 (including its sign error)
  Almkn = - 2.0 * II * omega * snt->Bin;

  // set integration boundaries
  IntToPi = 1;

  // Turn integrand scaling on or off.  If on, the radial integrand is rescaled
  // by a change of varriables according to:
  //
  //  int_0^pi dpsi f(psi)  -->  int_A^B dy exp(-psi/ScaleLambda) f(psi)
  //
  if(IntToPi)
    ScaleIntegrand = 1;
  else
    ScaleIntegrand = 0;

  // compute alpha_{lmkn} addapted from RRGW::alpha_func([stuff])
  //
  // Factor that appears in the fluxes on the horizon; arises from
  // the conversion of the Newman-Penrose tetrad to a Hartle-
  // Hawking tetrad (used to describe how the horizon grows due to
  // stuff falling in).  See Teukolsky and Press 1974.
  //
  const Real lamb = swsh->lambda;
  const Real w = omega;
  const Real r_plus = Kerr::rplus(a);
  const Real P = w - Real(m)*a/(2.0*r_plus);
  const Real aeps = sqrt(1. - a*a)/(4.*r_plus);
  const Real lp2 = lamb + 2.;
  const Real aw_m_m = a*w - ((Real)m);
  const Real tmp1 = lp2*lp2 - 4.*a*w*aw_m_m;
  const Real tmp2 = lamb*lamb - 36.*a*w*aw_m_m;
  const Real tmp3 = 48.*a*w*(2.*lamb + 3.)*(2.*a*w - ((Real)m));
  const Real tmp4 = 256.*pow(2.*r_plus, 5.)*P*(P*P + 4.*aeps*aeps)*
    (P*P + 16.*aeps*aeps)*w*w*w;
  const Real abssqrClm = tmp1*tmp2 + tmp3 + 144.*w*w*(1. - a*a);
  alpha = tmp4/abssqrClm;

  // enter this loop if we want to compute ZH
  if(DoZH){
    // determine rescaling rules for the Z^H_{lmkn} integrand
    DoH = 1;
    EstimateScaleLambda();
    Real A = 0.0;
    Real B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

    // compute Z^H_{lmkn}
    int PsiIntegral=1;
    if(ScaleIntegrand) {
      ZH = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, Min(1e-5,eps));
    } else {
      if(IntToPi)
        ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, Min(1e-5,eps));
      else
        ZH = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, Min(1e-5,eps));
    }
    ZH /= -Almkn;
    ZH *= 1.0 / (2.0 * M_PI * Gamma);
  } else {
   ZH = 0.0;
  }

  // enter this loop if we want to compute ZInf
  if(DoZI){
    // determine rescaling rules for the Z^inf_{lmkn} integrand
    DoH = 0;
    EstimateScaleLambda();
    Real A = 0.0;
    Real B = ScaleLambda*(exp(M_PI/ScaleLambda) - 1.0);

    // compute Z^{\infty}_{lmkn}
    int PsiIntegral=1;
    if(ScaleIntegrand) {
      ZInf = Integrate(&IEKR::OuterPsiIntegrand,A,B,PsiIntegral, Min(1e-5,eps));
    } else {
      if(IntToPi)
        ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,M_PI,PsiIntegral, Min(1e-5,eps));
      else
        ZInf = Integrate(&IEKR::OuterPsiIntegrand,0.0,2.0*M_PI,PsiIntegral, Min(1e-5,eps));
    }
    ZInf /= -Almkn;
    ZInf *= 1.0 / (2.0 * M_PI * Gamma);
  } else {
    ZInf = 0.0;
  }

}

// integrand for psi integral
Complex IEKR::OuterPsiIntegrand(Real psi_tmp){

  // simple quantities to be passed around
  if(ScaleIntegrand)
    psi = ScaleLambda * log(1.0 + psi_tmp/ScaleLambda);
  else
    psi = psi_tmp; 

  // this integrand is only be defined for 0 < psi < pi
  if( (psi < 0 || psi > M_PI) && IntToPi){
    if( psi < 0.0 && fabs(psi) < 1e-6)
      psi = 0.0;
    else if(psi > M_PI && fabs(psi-M_PI) < 1e-6)
      psi = M_PI;
    else{
      cout << "Psi  = " << psi << endl;
      Die("Psi out of range in IEKR::PsiIntegrand");
    }
  }

  // compute r from psi
  r = iekg->RofPsi(psi); 

  // K and Delta stuff
  Delta = Kerr::Delta(r,a);
  Real dDelta = Kerr::dr_Delta(r);
  Real K = (r*r + a*a)*omega - Real(m)*a;
  Real dK = 2.0*r*omega;
  KoverDelta = K/Delta;
  dKoverDelta =  (Delta*dK - K*dDelta) / (Delta*Delta);

  // compute RDot = dr/dlambda 
  Real RDotSquared = gkg.RFunc(r, a, 0.0, E, Lz, Q); // third argument unused
  RDot = ( RDotSquared > 0.0 ?  sqrt(RDotSquared) : 0.0);
  //
  // or use this
  //
  //Real RDotSquared = (E*(r*r + a*a) - a*Lz ) * (E*(r*r + a*a) - a*Lz );
  //RDotSquared -= Delta*(r*r + (Lz-a*E)*(Lz-a*E) + Q);
  //RDot = ( RDotSquared > 0.0 ?  sqrt(RDotSquared) : 0.0);

  // compute the homogeneous solution at new radius
  snt->CalcRXFields(r,DoH);
  //snt->SplineFields(r,DoH);
  if(DoH){
    R = snt->TeukRH;
    dR = snt->dr_TeukRH;
    ddR = snt->ddr_TeukRH;
  } else {
    R = snt->TeukRInf;
    dR = snt->dr_TeukRInf;
    ddR = snt->ddr_TeukRInf;
  }

  // do the chi integral
  Complex output = iekg->dwdpsi(psi);
  int PsiIntegral = 0;
  if(IntToPi)
    output *=  Integrate(&IEKR::InnerChiIntegrand,0.0,M_PI,PsiIntegral, Min(eps/10.0,1e-6));
  else
    output *=  Integrate(&IEKR::InnerChiIntegrand,0.0,2.0*M_PI,PsiIntegral, Min(eps/10.0,1e-6));


  // scale the integrand if necessary
  if(ScaleIntegrand) output *= exp(-psi/ScaleLambda);

  // exit
  return output;

}

// integrand for chi integral
Complex IEKR::InnerChiIntegrand(Real chi_tmp){

  // this integrand is only be defined for 0 < chi < pi
  if( (chi_tmp < 0 || chi_tmp > M_PI) && IntToPi) 
    Die("Chi out of range in IEKR::ChiIntegrand");

  // simple quantities to be passed around
  chi = chi_tmp;
  theta = iekg->Theta(chi);
  ct = cos(theta);
  st = sin(theta);
  cott = ct/st;
  csct = 1.0/st;
  Sigma = Kerr::Sigma(r,a,ct*ct);
  rho = -1.0 / (r - II*a*ct);
  rhobar = conj(rho);
  wR = iekg->wR(psi);
  wTheta = iekg->wTheta(chi);

  // lambda derivatives of t and theta
  TDot = gkg.TFunc(r, a, ct*ct, E, Lz, Q);
  Real ThetaDotSquared = gkg.ThetaFunc(0.0, a, ct*ct, E, Lz, Q); // first argument unused
  ThetaDot = (ThetaDotSquared > 0.0 ? sqrt(ThetaDotSquared) : 0.0);
  //
  // or use this
  //
  //TDot = iekg->TofPsi(psi) + iekg->TofChi(chi);
  //Real ThetaDotSquared = Q - cott*cott*Lz*Lz - a*a*ct*ct*(1.0 - E*E);
  //ThetaDot = (ThetaDotSquared > 0.0 ? sqrt(ThetaDotSquared) : 0.0);

  // spheroidal harmonic and its derivatives
  S = swsh->spheroid(ct);
  L2S = swsh->l2dagspheroid(ct);
  L1L2S = swsh->l1dagl2dagspheroid(ct);

  // the exclusively-chi part of the integral
  Complex output = iekg->dwdchi(chi);

  if(IntToPi){
    // sum the different branches of the integrand
    Complex branches = (0.0,0.0);
    Complex ExpFactor;
    for(SgnRDot = -1; SgnRDot <= 1; SgnRDot += 2){
      for(SgnThetaDot = -1; SgnThetaDot <= 1; SgnThetaDot += 2){
  
        // give derivatives proper signs
        RDot = SgnRDot * fabs(RDot);
        ThetaDot = SgnThetaDot * fabs(ThetaDot);
          
        // main part of the integral
        ExpFactor = exp( II * (SgnRDot*n*wR + SgnThetaDot*k*wTheta)  );
        branches += ExpFactor * Jlm();
      }
    } 
    // scale output by the branches
    output *= branches;
  } else {
    SgnRDot = (psi < M_PI ? 1 : -1);
    SgnThetaDot = (chi < M_PI ? 1 : -1);
    RDot = SgnRDot * fabs(RDot);
    ThetaDot = SgnThetaDot * fabs(ThetaDot);
    output *= Jlm() * exp( II * (n*wR + k*wTheta) );
  }

  // 3D plot test 
  /*
  if(DoH){
    cout.precision(16);
    Real yy = ScaleLambda*(exp(psi/ScaleLambda) - 1.0);
    if(ScaleIntegrand)
      cout << yy << "\t" << chi << "\t" << exp(-psi/ScaleLambda)*abs( (iekg->dwdpsi(psi))*output ) << "\n";
    else
      cout << psi << "\t" << chi << "\t" << abs( (iekg->dwdpsi(psi))*output ) << "\n";
  }
  */

  // exit
  return output;
}

// This is the raw integrand for the Z_{lmkn}'s.  It's used mainly for probing 
// general behavior (for example, for determining how to rescale the radial integral).
//
// Z_{lmnk} = (2\pi\Gamma)^{-1} \int_0^{2\pi} d\psi  \int_0^{2\pi} d\chi RawIntegrand(chi,psi)
Complex IEKR::RawIntegrand(Real chi_tmp, Real psi_tmp){

  // get psi and chi from input
  chi = chi_tmp;
  psi = psi_tmp;

  // this integrand is only be defined for 0 < (psi,chi) < pi
  if( (psi < 0 || psi > M_PI) && IntToPi){
    if( psi < 0.0 && fabs(psi) < 1e-12)
      psi = 0.0;
    else if(psi > M_PI && fabs(psi-M_PI) < 1e-12)
      psi = M_PI;
    else
      Die("Psi out of range in IEKR::PsiIntegrand");
  }
  if( (chi < 0 || chi > M_PI) && IntToPi)
    Die("Chi out of range in IEKR::ChiIntegrand");

  // begin to construct the output
  Complex output = iekg->dwdpsi(psi);
  output *= iekg->dwdchi(chi);

  // simple quantities to be passed around
  r = iekg->RofPsi(psi);
  theta = iekg->Theta(chi);
  ct = cos(theta);
  st = sin(theta);
  cott = ct/st;
  csct = 1.0/st;
  Sigma = Kerr::Sigma(r,a,ct*ct);
  rho = -1.0 / (r - II*a*ct);
  rhobar = conj(rho);
  wR = iekg->wR(psi);
  wTheta = iekg->wTheta(chi);

  // K and Delta stuff
  Delta = Kerr::Delta(r,a);
  Real dDelta = Kerr::dr_Delta(r);
  Real K = (r*r + a*a)*omega - Real(m)*a;
  Real dK = 2.0*r*omega;
  KoverDelta = K/Delta;
  dKoverDelta =  (Delta*dK - K*dDelta) / (Delta*Delta);

  // compute |dr/dlambda|
  Real RDotSquared = gkg.RFunc(r, a, 0.0, E, Lz, Q); // third argument unused
  RDot = ( RDotSquared > 0.0 ?  sqrt(RDotSquared) : 0.0);

  // compute |dtheta/dlambda|
  Real ThetaDotSquared = gkg.ThetaFunc(0.0, a, ct*ct, E, Lz, Q); // first argument unused
  ThetaDot = (ThetaDotSquared > 0.0 ? sqrt(ThetaDotSquared) : 0.0);

  // lambda derivative of t
  TDot = gkg.TFunc(r, a, ct*ct, E, Lz, Q);

  // compute the homogeneous solution at new radius
  snt->CalcRXFields(r,DoH);
  if(DoH){
    R = snt->TeukRH;
    dR = snt->dr_TeukRH;
    ddR = snt->ddr_TeukRH;
  } else {
    R = snt->TeukRInf;
    dR = snt->dr_TeukRInf;
    ddR = snt->ddr_TeukRInf;
  }
    
  // spheroidal harmonic and its derivatives
  S = swsh->spheroid(ct);
  L2S = swsh->l2dagspheroid(ct);
  L1L2S = swsh->l1dagl2dagspheroid(ct);

  if(IntToPi){
    // sum the different branches of the integrand
    Complex branches = (0.0,0.0);
    Complex ExpFactor;
    for(SgnRDot = -1; SgnRDot <= 1; SgnRDot += 2){
      for(SgnThetaDot = -1; SgnThetaDot <= 1; SgnThetaDot += 2){

        // give derivatives proper signs
        RDot = SgnRDot * fabs(RDot);
        ThetaDot = SgnThetaDot * fabs(ThetaDot);

        // main part of the integral
        ExpFactor = exp( II * (SgnRDot*n*wR + SgnThetaDot*k*wTheta)  );
        branches += ExpFactor * Jlm();
      }
    }
    // scale output by the branches
    output *= branches;
  } else {
    SgnRDot = (psi < M_PI ? 1 : -1);
    SgnThetaDot = (chi < M_PI ? 1 : -1);
    RDot = SgnRDot * fabs(RDot);
    ThetaDot = SgnThetaDot * fabs(ThetaDot);
    output *= Jlm() * exp( II * (n*wR + k*wTheta) );
  }

  // exit
  return output;

}

// The function EstimateScaleLambda sets a value for ScaleLambda
// (it can also set ScaleIntegrand = 0 if it doesn't look like there's any
// point to rescaling).
void IEKR::EstimateScaleLambda(){

  // if rescaling has been turned off, don't waste time in here.
  // Just set ScaleLambda to some arbitrary big-ish number (which won't be used
  // anyway) and return.
  if(ScaleIntegrand != 1){
    ScaleLambda = 100.0*M_PI;
    return;
  }
  
  // get average |integrand| at psi = pi
  Real Ipi = 0;
  Real dummypsi = M_PI;
  for(Real dummychi = 0.0; dummychi < M_PI; dummychi += M_PI/15.0){
    Ipi += abs(RawIntegrand(dummychi, dummypsi)) / 15.0;
  }
  if(Ipi == 0){ // something went wrong, shut off rescaling
    cout << "%WARNING IEKR::EstimateScaleLambda: Ipi = 0, rescaling is shut off.\n";
    ScaleLambda = 100.0*M_PI;
    ScaleIntegrand = 0;
    return;
  }

  // get average |integrand|  psi = 0
  Real I0 = 0;
  dummypsi = 0;
  for(Real dummychi = 0.0; dummychi < M_PI; dummychi += M_PI/15.0){
    I0 += abs(RawIntegrand(dummychi, dummypsi)) / 15.0;
  }
  if(I0 == 0){ // something went wrong, shut off rescaling
    cout << "%WARNING IEKR::EstimateScaleLambda: I0 = 0, rescaling is shut off.\n";
    ScaleLambda = 100.0*M_PI;
    ScaleIntegrand = 0;
    return;
  }

  // Set value for the estimated rescaling parameter.
  ScaleLambda = M_PI / log(Ipi/I0);

  // We only rescale if things Ipi and I0 differ by more than a factor of 10
  if( exp(-M_PI/ScaleLambda) < 10 && exp(-M_PI/ScaleLambda) > 0.1) {
    ScaleLambda = 100.0*M_PI;
    ScaleIntegrand = 0;
  }

  // if we are damping by 10 orders of magnitude or more, it's warning time.
  if( exp(-M_PI/ScaleLambda) > 1e10 || exp(-M_PI/ScaleLambda) < 1e-10)
    cout << "%WARNING IEKR::EstimateScaleLambda: Extreme damping requested!\n";

  // all done
  return;

}

Complex IEKR::Jlm()
{

  // compute the exponental and dt/d\lambda prefactor
  Complex output = TDot;
  if(IntToPi) {
    output *= exp( II* 
                   ( omega*iekg->Delta_t(SgnRDot,SgnThetaDot,wR,wTheta) 
                       - m*iekg->Delta_phi(SgnRDot,SgnThetaDot,wR,wTheta) )  
                 );
  } else {
    output *= exp( II*
                   ( omega*iekg->Delta_t(wR,wTheta)
                       - m*iekg->Delta_phi(wR,wTheta) )
                 );
  }


  // compute I_{lm}(r,theta,omega)
  output *= Ilm();

  // exit
  return output;
}
                                                                                                              
// I_{lm}(r,theta,omega)
Complex IEKR::Ilm(){

  // terms propto R
  Complex output = ( A_n_n_0() + A_n_mbar_0() + A_mbar_mbar_0() ) * R;

  // terms propto d/dr  R
  output -= ( A_n_mbar_1() + A_mbar_mbar_1() ) * dR;

  // term propto d/dr d/dr R
  output += A_mbar_mbar_2() * ddR;

  // exit
  return output;
}

// the A_{abc}
inline Complex IEKR::A_n_n_0(){
  Complex output = L1L2S + 2.0*II*a*rho*st*L2S;
  output *= - 2.0 / ( rho*rho*rho * rhobar * Delta*Delta);
  output *= C_n_n();
  return output;
}
inline Complex IEKR::A_n_mbar_0(){
  Complex output = (II*KoverDelta - rho - rhobar) * L2S;
  output += (II*KoverDelta + rho + rhobar) * II*a*st*S*(rho-rhobar);
  output *= -2.0*sqrt(2.0) / (rho*rho*rho * Delta);
  output *= C_n_mbar(); 
  return output;
}
inline Complex IEKR::A_mbar_mbar_0(){
  Complex output = KoverDelta*KoverDelta;
  output += 2.0*II*rho*KoverDelta;
  output += II*dKoverDelta;
  output *= S*rhobar / (rho*rho*rho);
  output *= C_mbar_mbar(); // this factor is missing from PRD 61 084004 
  return output;
}
inline Complex IEKR::A_n_mbar_1(){
  //Complex output = L2S + II*a*rho*st*S*(rho - rhobar); // mistake in PRD 61 084004
  Complex output = L2S + II*a*st*S*(rho - rhobar);  // fix
  output *= -2.0*sqrt(2.0) / (rho*rho*rho * Delta);
  output *= C_n_mbar();
  return output;
}
inline Complex IEKR::A_mbar_mbar_1(){
  Complex output = 2.0*S*rhobar*(rho - II*KoverDelta) / (rho*rho*rho);
  output *= C_mbar_mbar();
  return output;
}
inline Complex IEKR::A_mbar_mbar_2(){
  Complex output = -S*rhobar / (rho*rho*rho);
  output *= C_mbar_mbar();
  return output;
}

// the C_{ab}'s
inline Real IEKR::C_n_n(){ 
  Real output = E*(r*r + a*a) - a*Lz + RDot;
  output *= output;
  output /= 4.0*TDot*Sigma*Sigma;
  return output;
}
inline Complex IEKR::C_mbar_mbar(){
  Complex output = II*(a*E - Lz*csct*csct)*st + ThetaDot;
  output *= output;
  output *= rho*rho / (2.0*TDot);
  return output;
}
inline Complex IEKR::C_n_mbar(){
  Complex output = II*(a*E - Lz*csct*csct)*st + ThetaDot;
  output *= E*(r*r + a*a) - a*Lz + RDot;
  output *= rho / (2.0*sqrt(2.0)*TDot*Sigma);
  return output;
}

// integration wrapper function
//
// There are several integration options available. Uncomment only one of them.
Complex IEKR::Integrate(Complex (IEKR::*igrnd)(Real), Real start, Real finish, int PsiIntegral, Real eps_int){

  // Clenshaw-Curtis: Uses analytical expression of integral for an integrand which
  //                  has been expressed as a Chebyshev series.  Uses rational function
  //                  interpolation to extrapolate to an infinite expansion (after a minimum
  //                  of 2 iterations).
  //
  // starts at N = ncofs terms, and doubles N on each iteration.
  // 2.^(4:15) = 16     32     64    128    256    512   1024   2048   4096   8192  16384  32768
  if(PsiIntegral == 1){
    int ncofs = 16;
    int itmax = 8; 
    return ClenshawCurtisInt(igrnd, start, finish, ncofs, itmax, eps_int);
  } else {
    int ncofs = 16;
    int itmax = 7;
    return ClenshawCurtisInt(igrnd, start, finish, ncofs, itmax, eps_int);
  }

  // Cash-Carp: An adaptive stepsize integration routine which uses (embedded) 
  //            4th and 5th order Runge-Kutta approximations to determine stepsize.
  //
  // Set minimum step size for adaptive step integrations:
  // guess 1000 steps per perriod
  /*
  Real h_guess;
  if(PsiIntegral == 1)
    h_guess = ( n==0 ? 2.0*M_PI/1000.0 : 2.0*M_PI / (1000.0*n) );
  else
    h_guess = ( k==0 ? 2.0*M_PI/1000.0 : 2.0*M_PI / (1000.0*k) );
  return CashCarpInt(igrnd, start, finish, h_guess, eps_int);
  */

  //
  // Romberg method:  Just the trapazoid rule, but uses polynomial extrapolation
  //                  to guess result for zero stepsize. 
  //
  // Set minimum number of iterations for Romberg integrations.
  // Try for 20 steps per perriod.
  /*
  if(PsiIntegral == 1)
    MinIterate = ( n==0 ? 5 : int(log(10.0*abs(n))/log(2.0) + 1) );
  else
    MinIterate = ( k==0 ? 5 : int(log(10.0*abs(k))/log(2.0) + 1) );
  return RombergInt(igrnd, start, finish, MinIterate, eps_int);
  */

}

// 
// Romberg integration routines
//
//  These are modified versions of the numerical recipes routines
//  qromb, trpzd, and polint.  There are three modifications:
//
//  1) don't use the numerical recipes vector class
//  2) eliminate the sneaky static varriable in trpzd
//  3) accomadate for complex functions
// 
Complex IEKR::RombergInt(Complex (IEKR::*func)(Real), Real start, Real finish, int MinIterate, Real eps_int)
{
        int JMAX=100, JMAXP=JMAX+1, K=MinIterate;
        Complex ss,dss;
        Complex *s, *s_t;
        Real *h, *h_t;
        s = Complexvector(0,JMAX-1);
        h = Realvector(0,JMAXP-1);
        s_t = Complexvector(0,K-1);
        h_t = Realvector(0,K-1);
        int i,j;

        h[0]=1.0;
        for (j=1;j<=JMAX;j++) {
                if(j==1){
                  Ctrpzd(func,start,finish,j,0.0,&s[0]);
                } else {
                  Ctrpzd(func,start,finish,j,s[j-2],&s[j-1]);
                }
                if (j >= K) {
                        for (i=0;i<K;i++) {
                                h_t[i]=h[j-K+i];
                                s_t[i]=s[j-K+i];
                        }
                        Cpolint(h_t,s_t,0.0,ss,dss,K);
                        //if (fabs(dss) <= eps_int*fabs(ss)) return ss;
                        //modified from above to account for complex
                        //numbers and 0 result
                        if (fabs(dss.real()) <= eps_int*fabs(ss.real()) &&
                           fabs(dss.imag()) <= eps_int*fabs(ss.imag()) ||
                           fabs(ss.real()) < eps_int && 
                           fabs(ss.imag()) < eps_int) {
                        // scott uses this
                        //if (fabs(dss.real()) <= eps_int*fabs(ss.real()) &&
                        //   fabs(dss.imag()) <= eps_int*fabs(ss.imag()) ){
                          free_Complexvector(s,0,JMAX-1);
                          free_Realvector(h,0,JMAXP-1);
                          free_Complexvector(s_t,0,K-1);
                          free_Realvector(h_t,0,K-1);
                          return ss;
                        }
                }
                h[j]=0.25*h[j-1];
        }
        Die("Too many steps in routine IEKR::Cqromb");
        return 0.0;
}
void IEKR::Ctrpzd(Complex (IEKR::*func)(Real), Real start, Real finish,
     int nn, Complex last, Complex *next)
{
        Real x,tnm,del;
        Complex sum;
        int it,j;

        if (nn == 1) {
                *next=0.5*(finish-start)*( (this->*func)(start) + (this->*func)(finish) );
                return;
        } else {
                for (it=1,j=1;j<nn-1;j++) it <<= 1;
                tnm=it;
                del=(finish-start)/tnm;
                x=start+0.5*del;
                for (sum=0.0,j=0;j<it;j++,x+=del) sum += (this->*func)(x);
                *next=0.5*(last+(finish-start)*sum/tnm);
                return;
        }
}
void IEKR::Cpolint(Real *xa, Complex *ya, Real x, Complex &y, Complex &dy, int nn)
{
        int i,mm,ns=0;
        Real dif,dift,ho,hp;
        Complex *c, *d, den, w;
        c = Complexvector(0,nn-1);
        d = Complexvector(0,nn-1);
        dif=fabs(x-xa[0]);
        for (i=0;i<nn;i++) {
                if ((dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        y=ya[ns--];
        for (mm=1;mm<nn;mm++) {
                for (i=0;i<nn-mm;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+mm]-x;
                        w=c[i+1]-d[i];
                        if (ho-hp == 0.0) Die("Divide by zero error in routine Cpolint");
                        den=w/(ho - hp);
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                y += (  dy=(2*(ns+1) < (nn-mm) ? c[ns+1] : d[ns--])  );
        }
        free_Complexvector(c,0,nn-1);
        free_Complexvector(d,0,nn-1);
}

//
// Adaptive stepsize integration routines
//
// based on odeint, rkqs, and rkck, in Numerical Recipes
//
Complex IEKR::CashCarpInt(Complex (IEKR::*derivs)(Real), Real x1, Real x2, Real h_guess, Real eps_int)
{
        const int MAXSTP=1000000;
        const Real TINY=1.0e-30;
        Real hmin = 0;
        int i,nstp;
        Real x,hnext,hdid,h;
        Real h1 = h_guess;

        Real yscal;
        Complex y;
        Complex dydx;
        x=x1;
        h=SIGN(h1,x2-x1);
        y = 0.0;
        for (nstp=0;nstp<MAXSTP;nstp++) {
                dydx=(this->*derivs)(x);
                yscal=abs(y)+abs(dydx*h)+TINY;
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                Crkqs(y,dydx,x,h,eps_int,yscal,hdid,hnext,derivs);
                if ((x-x2)*(x2-x1) >= 0.0) {
                        return y;
                }
                if (fabs(hnext) <= hmin) Die("Step size too small in IEKR::Codeint");
                h=hnext;
        }
        Die("Too many steps in routine IEKR::Codeint");
        return 0.0; // never make it here. this just quiets a compile warning
}
void IEKR::Crkqs(Complex &y, Complex dydx, Real &x, Real htry,
        Real eps_int, Real yscal, Real &hdid, Real &hnext,
        Complex (IEKR::*derivs)(Real) )
{
        const Real SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
        int i;
        Real errmax,h,htemp,xnew;

        h=htry;
        Complex yerr;
        Complex ytemp;
        for (;;) {
                Crkck(y,dydx,x,h,ytemp,yerr,derivs);
                errmax=0.0;
                errmax=Max(errmax,abs(yerr/yscal));
                errmax /= eps_int;
                if (errmax <= 1.0) break;
                htemp=SAFETY*h*pow(errmax,PSHRNK);
                h=(h >= 0.0 ? Max(htemp,0.1*h) : Min(htemp,0.1*h));
                xnew=x+h;
                if (xnew == x) Die("stepsize underflow in IEKR::Crkqs");
        }
        if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
        else hnext=5.0*h;
        x += (hdid=h);
        y=ytemp;
}
void IEKR::Crkck(Complex &y, Complex &dydx, Real x,
        Real h, Complex &yout, Complex &yerr,
        Complex (IEKR::*derivs)(Real) )
{
        static const Real a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
                b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
                b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
                b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
                b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
                c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
                dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
                dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
        int i;

        Complex ak2, ak3, ak4, ak5, ak6;
        ak2=(this->*derivs)(x+a2*h);
        ak3=(this->*derivs)(x+a3*h);
        ak4=(this->*derivs)(x+a4*h);
        ak5=(this->*derivs)(x+a5*h);
        ak6=(this->*derivs)(x+a6*h);
        yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6);
        yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);
}

//
// Clenshaw-Curtis Integration routines (adapted from CKR.cc so 
// as to be used on arbitrary integrals)
//
Complex IEKR::ClenshawCurtisInt(Complex (IEKR::*igrnd)(Real), Real x1, Real x2, int ncofs, int itmax, Real eps_int){
  int isradial = (itmax == 8 ? 1 : 0);
  int isangular = (itmax == 11 ? 1 : 0);
  int iter;
  int i;
  Complex sum;
  Real *oldy, *newy;
  Real pio2 = 0.5*M_PI;
  Complex *oldcofs, *newcofs;
  Complex maxcof;
  //
  // Storage and utility variables for the rational extrapolation.
  //
  int K = 3;
  Complex Sum, Sumerr;
  Complex *sumest;    sumest = Complexvector(1, itmax);
  Real *oneovern;      oneovern = Realvector(1, itmax);
  for (iter = 1; iter <= itmax; iter++) {
    if (iter == 1) {
      newy = Realvector(1, ncofs + 1);
      oldy = Realvector(1, ncofs + 1);
      newcofs = Complexvector(1, ncofs + 1);
      oldcofs = Complexvector(1, ncofs + 1);
      for (i = 1; i <= ncofs + 1; i++) {
        //
        // Evaluating the integrand at the proper points in prep for
        // the Chebyshev expansion.  See NR (2nd ed), Eqs (5.9.4) and
        // (5.8.10) (with a = 0, b = pi). 
        //
        newy[i] = cos(M_PI*(i - 1)/((Real)ncofs));
        //newcofs[i] = (this->*igrnd)(newy[i]*pio2 + pio2); // only for int_0^pi
        newcofs[i] = (this->*igrnd)(0.5*(x2-x1)*newy[i] + 0.5*(x2+x1));  // for generic boundaries
        //
        // Copy data into "old" arrays to be used in the next
        // iteration.
        //
        oldy[i] = newy[i];
        oldcofs[i] = newcofs[i];
      }
      //
      // The cosine Fourier transform.  Note, I have modified the NR
      // FFT routines very slightly to handle complex functions.
      // Since there is no crosstalk between the real and the
      // imaginary parts, this is equivalent to FFT'ing a pair of real
      // functions.  I found experimentally that doing it this way is
      // quite a bit faster than storing the real and imaginary parts
      // separately and performing two FFTs.
      //
      cosft1(newcofs, ncofs);
      //
      // Clenshaw-Curtis quadrature --- NR ed 2, Eq (5.9.3).
      maxcof = newcofs[1];
      sum = 0.5*newcofs[1];
      for (i = 1; i <= ncofs/2; i++) {
        sum -= newcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
        //
        if (abs(newcofs[2*i + 1]) > abs(maxcof))
          maxcof = newcofs[2*i + 1];
      }
    } else {
      //
      // This bracketed section is identical to the above section,
      // except that we carefully reuse data that was calculated on
      // the previous iteration whenever possible.
      //
      ncofs *= 2;
      newy = Realvector(1, ncofs + 1);
      newcofs = Complexvector(1, ncofs + 1);
      //
      // Points we can recycle from the previous iteration
      for (i = 1; i <= ncofs + 1; i += 2) {
        newy[i] = oldy[(i - 1)/2 + 1];
        newcofs[i] = oldcofs[(i - 1)/2 + 1];
      }
      //
      // Points we cannot recycle
      for (i = 2; i <= ncofs; i += 2) {
        newy[i] = cos(M_PI*(i - 1)/((Real)ncofs));
        //newcofs[i] = (this->*igrnd)(newy[i]*pio2 + pio2); // only for int_0^pi
        newcofs[i] = (this->*igrnd)(0.5*(x2-x1)*newy[i] + 0.5*(x2+x1)); // for generic boundaries
      }
      //
      // Remake the "old" arrays
      free_Realvector(oldy, 1, ncofs/2 + 1);
      oldy = Realvector(1, ncofs + 1);
      free_Complexvector(oldcofs, 1, ncofs/2 + 1);
      oldcofs = Complexvector(1, ncofs + 1);
      for (i = 1; i <= ncofs + 1; i++) {
        oldy[i] = newy[i];
        oldcofs[i] = newcofs[i];
      }
      cosft1(newcofs, ncofs);
      //
      maxcof = newcofs[1];
      sum = 0.5*newcofs[1];
      for (i = 1; i <= ncofs/2; i++) {
        sum -= newcofs[2*i + 1]/((Real)((2*i + 1)*(2*i - 1)));
        if (abs(newcofs[2*i + 1]) > abs(maxcof))
          maxcof = newcofs[2*i + 1];
      }
    }
    //sum *= 2.*M_PI/((Real)ncofs); // only for int_0^pi
    //maxcof *= 2.*M_PI/((Real)ncofs); // only for int_0^pi
    sum *= 2.*(x2-x1)/((Real)ncofs); // for generic boundaries
    maxcof *= 2.*(x2-x1)/((Real)ncofs); // for generic boundaries
    free_Realvector(newy, 1, ncofs + 1);
    free_Complexvector(newcofs, 1, ncofs + 1);
    //
    // This threshold of 1e-8 seems to be somewhat robust for the IE code.  
    // The circular code has been able to get away with a much smaller 
    // threshold of around 1e-11.
    if ( abs(sum) < (1.e-7)*abs(maxcof) )  {
      // The sum is either zero, or so close to zero we're just getting roundoff crap.
      // Set it to zero and break.
      Sum *= 0.0;
      break;
    }
    sumest[iter] = sum;

    // In order to better represent the way we actually chose the number of 
    // coefficients, we use a log-scale x-axis when extrapolating to an infinite number
    // of coefficients.
    oneovern[iter] = fabs(  1./log2(Real(ncofs))   );
    
    // fit a rational function to estimate(1/N).  Converge when 
    // the error estimate on an extrapolation to estimate(0) is 
    // sufficiently small
    if (iter >= K) {
      ratint(&oneovern[iter - K], &sumest[iter - K], K, 0., &Sum, &Sumerr);
      if (abs(Sumerr) <= eps_int*abs(Sum) ) break;
    }
    // We'll also be happy to converge when fractional changes 
    // in consecutive estimates are sufficiently small
    if (iter>1){
      Real FracChange = abs(sumest[iter-1] - sumest[iter]) / abs(sumest[iter]);
      if( FracChange <= eps_int  ) {
        Sum = sum;
        break;
      }
    }

    // If we know this integral to an accuracy such that the uncertainty
    // in the energy flux will be negligable, then we stop iterating.
    //
    if(isradial && iter >= K && abs(Sumerr)/abs(Sum) < 1e-1){
      Complex ZEstimate = -Sum/(2.0*M_PI*Gamma*Almkn);
      Real EfluxEst = abs(ZEstimate)*abs(ZEstimate) / (4.0 * M_PI * omega * omega);
      Real epsc = abs(Sumerr) /abs(Sum);
      // we put in an extra saftey factor
      Real saftey = 1.e-1;
      if(  (DoH == 1 && 2.0*epsc*EfluxEst             < saftey*eps_int*MaxEI) || 
           (DoH == 0 && 2.0*epsc*fabs(alpha*EfluxEst) < saftey*eps_int*MaxEH)  ) break;
    }

  } // end of main iteration loop of the integrator



  // We might need a warning if this integral failed to converge.
  if(isradial && iter >= itmax){ 
    // A radial integral didn't converge.  
    //
    // Estimate whether or not it was significant.  If it was, write a warning.
    Complex ZEstimate = -Sum/(2.0*M_PI*Gamma*Almkn);
    Real EfluxEst = abs(ZEstimate)*abs(ZEstimate) / (4.0 * M_PI * omega * omega);
    Real epsc = abs(Sumerr) /abs(Sum);
    Real saftey = 1.e-1;
    int Significant = 0;
    if(  (DoH == 1 && 2.0*epsc*EfluxEst             > saftey*eps_int*MaxEI) ||
         (DoH == 0 && 2.0*epsc*fabs(alpha*EfluxEst) > saftey*eps_int*MaxEH)  ) Significant = 1;
    if(Significant){
      cout.precision(16);
      cout << "% ---------------------------------------------------------\n"
           << "% WARNING a significant RADIAL INTEGRAL failed to converge!\n"
           << "% More than " << ncofs << " terms are needed.\n"
           << "% abs(sum)/abs(maxcof) = " << abs(sum)/abs(maxcof) << endl;
      if(DoH == 1) {
        cout << "% EI    = " << EfluxEst << endl;
        cout << "% MaxEI = " << MaxEI << endl;
      } else {
        cout << "% |EH|    = " << fabs(alpha*EfluxEst) << endl;
        cout << "% |MaxEH| = " << MaxEH << endl;
      }
      cout << "% eps_int = " << eps_int << endl;
      cout << "% Returning a zero result" << endl
           << "% ---------------------------------------------------------\n";
    }
    Sum *= 0.0;
  } else if(isangular && iter >= itmax){
      // An angular integral didn't converge
      // 
      // We don't have an estimate of signifigance for angular integrals.
      cout.precision(16);
      cout << "% ----------------------------------------------------------\n"
           << "% WARNING an ANGULAR INTEGRAL failed to converge!\n"
           << "% More than " << ncofs << " terms are needed.\n"
           << "% abs(sum)/abs(maxcof) = " << abs(sum)/abs(maxcof) << endl
           << "% eps_int = " << eps_int << endl;
      cout << "% Returning a zero result" << endl;
      cout << "% ----------------------------------------------------------\n";
      Sum *= 0;
  } // end of warning block

  // test (shows integral result as a function of iteration number)
  /*
  for(int ind = 1; ind <= iter; ind++){
    if(ind <= itmax) cout << ind << "\t" << sumest[ind] << endl;
  }
  cout << "---\n";
  */

  free_Realvector(oldy, 1, ncofs + 1);
  free_Complexvector(oldcofs, 1, ncofs + 1);
  free_Complexvector(sumest, 1, itmax);
  free_Realvector(oneovern, 1, itmax);
  return Sum;
}
void IEKR::cosft1(Complex y[], int nn)
{
  int j, n2;
  Complex sum, y1, y2;
  Real theta, wi=0.0, wpi, wpr, wr=1.0, wtemp;
                                                                                                    
  theta = M_PI/nn;
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  sum = 0.5*(y[1] - y[nn + 1]);
  y[1] = 0.5*(y[1] + y[nn + 1]);
  n2 = nn + 2;
  for (j = 2; j <= (nn>>1); j++) {
    wr = (wtemp = wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
    y1 = 0.5*(y[j] + y[n2 - j]);
    y2 = (y[j] - y[n2 - j]);
    y[j] = y1 - wi*y2;
    y[n2 - j] = y1 + wi*y2;
    sum += wr*y2;
  }
  realft(y, nn, 1);
  y[nn + 1] = y[2];
  y[2] = sum;
  for (j = 4; j <= nn; j += 2) {
    sum += y[j];
    y[j] = sum;
  }
}
void IEKR::realft(Complex data[], unsigned long nn, int isign)
{
  unsigned long i, i1, i2, i3, i4, np3;
  Real c1 = 0.5, c2;
  Complex h1r, h1i, h2r, h2i;
  Real wr, wi, wpr, wpi, wtemp, theta;
                                                                                                    
  theta = 3.141592653589793/(Real) (nn>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, nn>>1, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = nn+3;
  for (i = 2; i <= (nn>>2); i++) {
    i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
    h1r = c1*(data[i1] + data[i3]);
    h1i = c1*(data[i2] - data[i4]);
    h2r = -c2*(data[i2] + data[i4]);
    h2i = c2*(data[i1] - data[i3]);
    data[i1] = h1r + wr*h2r - wi*h2i;
    data[i2] = h1i + wr*h2i + wi*h2r;
    data[i3] = h1r - wr*h2r + wi*h2i;
    data[i4] = -h1i + wr*h2i + wi*h2r;
    wr = (wtemp = wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
  }
  if (isign == 1) {
    data[1] = (h1r = data[1]) + data[2];
    data[2] = h1r - data[2];
  } else {
    data[1] = c1*((h1r = data[1]) + data[2]);
    data[2] = c1*(h1r - data[2]);
    four1(data, nn>>1, -1);
  }
}
void IEKR::four1(Complex data[], unsigned long nn, int isign)
{
  unsigned long n, mmax, m, j, istep, i;
  Real wtemp, wr, wpr, wpi, wi, theta;
  Complex tempr, tempi;
                                                                                                    
  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i+= 2) {
    if (j > i) {
      tempr = data[j];
      data[j] = data[i];
      data[i] = tempr;
      //
      tempr = data[j + 1];
      data[j + 1] = data[i + 1];
      data[i + 1] = tempr;
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(2.*M_PI/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m+= 2) {
      for (i = m; i <= n; i+= istep) {
        j = i + mmax;
        tempr = wr*data[j] - wi*data[j + 1];
        tempi = wr*data[j + 1] + wi*data[j];
        data[j] = data[i] - tempr;
        data[j+1] = data[i + 1] - tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}
void IEKR::ratint(const Real xa[], const Complex ya[], const int n,
                 const Real x, Complex *y, Complex *dy)
{
  int m, i, ns = 1;
  const Real TINY = 1.0e-25;
  Real hh, h;
  Complex w, t, dd, *c, *d;
                                                                                                    
  c = Complexvector(1, n);
  d = Complexvector(1, n);
  hh = fabs(x - xa[1]);
  for (i = 1; i <= n; i++) {
    h = fabs(x - xa[i]);
    if (h == 0.0) {
      *y = ya[i];
      *dy = 0.0;
      free_Complexvector(d, 1, n);
      free_Complexvector(c, 1, n);
      return;
    } else if (h < hh) {
      ns = i;
      hh = h;
    }
    c[i] = ya[i];
    d[i] = ya[i] + TINY;
  }
  *y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      w = c[i+1] - d[i];
      h = xa[i+m] - x;
      t = (xa[i] - x)*d[i]/h;
      dd = t - c[i+1] + TINY;
      if (dd.real() == 0.0 && dd.imag() == 0.0)
        Die("Error in routine IEKR::ratint");
      dd = w/dd;
      d[i] = c[i+1]*dd + TINY;
      c[i] = t*dd;
    }
    *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_Complexvector(d, 1, n);
  free_Complexvector(c, 1, n);
  return;
}


// code from NR for error function: erf(x) = (2/pi) int_0^x dy exp(-y^2)
Real IEKR::erff(Real x)
{
        return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}
Real IEKR::gammp(Real aa, Real x)
{
        Real gamser,gammcf,gln;
                                                                                                                                                             
        if (x < 0.0 || aa <= 0.0)
                Die("IEKR::gammp: Invalid arguments");
        if (x < aa+1.0) {
                gser(gamser,aa,x,gln);
                return gamser;
        } else {
                gcf(gammcf,aa,x,gln);
                return 1.0-gammcf;
        }
}
void IEKR::gser(Real &gamser, Real aa, Real x, Real &gln)
{
        const int ITMAX=100;
        Real epsilon = 1e-11;
        int nn;
        Real sum,del,ap;
                                                                                                                                                             
        gln=gammln(aa);
        if (x <= 0.0) {
                if (x < 0.0) Die("IEKR::gser: x less than 0");
                gamser=0.0;
                return;
        } else {
                ap=aa;
                del=sum=1.0/aa;
                for (nn=0;nn<ITMAX;nn++) {
                        ++ap;
                        del *= x/ap;
                        sum += del;
                        if (fabs(del) < fabs(sum)*epsilon) {
                                gamser=sum*exp(-x+aa*log(x)-gln);
                                return;
                        }
                }
                Die("IEKR::gser: aa too large, ITMAX too small");
                return;
        }
}
void IEKR::gcf(Real &gammcf, Real aa, Real x, Real &gln)
{
        const int ITMAX=100;
        Real epsilon = 1e-11;
        Real FPMIN=1e-15;
        int i;
        Real an,b,c,d,del,h;
                                                                                                                                                             
        gln=gammln(aa);
        b=x+1.0-aa;
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=ITMAX;i++) {
                an = -i*(i-aa);
                b += 2.0;
                d=an*d+b;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=b+an/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) <= epsilon) break;
        }
        if (i > ITMAX) Die("IEKR::gcf: aa too large, ITMAX too small");
        gammcf=exp(-x+aa*log(x)-gln)*h;
}
Real IEKR::gammln(Real xx)
{
        int j;
        Real x,y,tmp,ser;
        static const Real cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
                -0.5395239384953e-5};
                                                                                                                                           
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<6;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}

