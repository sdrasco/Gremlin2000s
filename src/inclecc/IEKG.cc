//
// Methods for inclined eccentric Kerr geodesics. 
// Taken from from CKG.cc
//
// [1] W. Schmidt Class. Quantum Grav. 19 2743 (2002)
// [2] S. A. Hughes, Phys. Rev. D 61 084004 (2000)
// [3] S. Drasco and S. A. Hughes, astro-ph/0308479 
//
// WARNING: the notation conventions in these references [1-3] are NOT the 
// same. For example, z_[2,3] = ( z_[1] )^2.  This code has a mix of the 
// different conventions.
// 
// Steve Drasco
//
#include "IEKG.h"
#include "NRIEKG.h"
#include "NRUtil.h"

//  constructor
IEKG::IEKG(const int ProGrade, const Real spin, const Real eccentricity,
	   const Real SemiLatusRectum, const Real ThetaMinus, const Real tol) :  
           prograde(ProGrade),   a(spin), e(eccentricity),
           p(SemiLatusRectum), theta_(ThetaMinus), eps(tol)
{

  // This code is for non-circular bound orbits. It can treat small 
  // eccentricity, but cannot use an exactly zero eccentricity.
  if (e == 0)
    Die("Zero eccentricity?  Use the circular code, punky.");

  // define z_ as in [1] -- NOT AS IN [2,3]
  z_ = cos(theta_);
  
  // define r1 (rmin) and r2 (rmax)
  rmin = p / (1.0 + e);
  rmax = p / (1.0 - e);
  Real r1 = rmin;
  Real r2 = rmax;
  
  // define Schmidt's functions f(r) etc. for r1
  Real Delta1 = r1*r1 - 2.0*r1 + a*a;
  Real f1 = r1*r1*r1*r1 + a*a*( r1*(r1+2.0) + z_*z_*Delta1 );
  Real g1 = 2.0*a*r1;
  Real h1 = r1*(r1 - 2.0) + z_*z_*Delta1/(1.0-z_*z_);
  Real d1 = (r1*r1 + a*a*z_*z_)*Delta1;
  
  // define Schmidt's functions f(r) etc. for r2
  Real Delta2 = r2*r2 - 2.0*r2 + a*a;
  Real f2 = r2*r2*r2*r2 + a*a*( r2*(r2+2.0) + z_*z_*Delta2 );
  Real g2 = 2.0*a*r2;
  Real h2 = r2*(r2 - 2.0) + z_*z_*Delta2/(1.0-z_*z_);
  Real d2 = (r2*r2 + a*a*z_*z_)*Delta2;
  
  // compute determinents
  Real kappa =          d1*h2 - d2*h1;
  Real epsilon =        d1*g2 - d2*g1;
  Real rho =            f1*h2 - f2*h1;
  Real eta =            f1*g2 - f2*g1;
  Real sigma =          g1*h2 - g2*h1;

  // compute energy 
  E = kappa*rho + 2.0*epsilon*sigma;
  E += pow(-1.0,double(prograde))*2.0*
    sqrt(sigma*(sigma*epsilon*epsilon + rho*epsilon*kappa
		- eta*kappa*kappa) );
  E /= rho*rho + 4.0*eta*sigma;
  E = sqrt(E);
  
  // compute Lz 
  //
  // for non-zero a, we can make an equation which is 
  // linear in Lz from Schmidt's Eq. (B.17), by scaling
  // the two equations and subtracting them so as to 
  // knock out the Lz^2 term to get  
  //
  //  Lz = (E*E*rho - kappa) / (2.0*E*sigma).
  //
  // This trick doesn't work in the a=0 limit, since 
  // then the "linear" equation doesn't involve Lz.
  // So, to allow for a=0, we use this ugly block
  // to get Lz...
  //
  if (a > 0.001){
    Lz = (E*E*rho - kappa) / (2.0*E*sigma);
  } else {
    Real BB = 2.0*g1*E / h1;
    Real CC = (d1 - f1*E*E) / h1;
    if(prograde){
      Lz = -0.5*BB + 0.5*sqrt(BB*BB - 4.0*CC);
    } else {
      Lz = -0.5*BB - 0.5*sqrt(BB*BB - 4.0*CC);
    }
  }

  // now compute Q
  Q = z_*z_*(  a*a*(1.0-E*E) + Lz*Lz/(1.0-z_*z_)  );
  
  // compute cos(iota) as in [2]
  CosIota = Lz/sqrt(Lz*Lz + Q);
  
  // some functions in [1] needed to compute frequencies
  Real A = a*a*(1.0 - E*E);
  Real B = - Lz*Lz - Q - a*a*(1.0 - E*E);
  zPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*A);
  zMinusSquared = z_*z_;
  zmzp = zMinusSquared / zPlusSquared;
  beta = a*sqrt(1.0 - E*E);
  Real betaOvera = sqrt(1.0 - E*E);
  aazPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*(1.0 - E*E));
  betaRootzPlusSquared = betaOvera*sqrt(aazPlusSquared);
  
  // some integrals in [1] needed to compute frequencies
  // Note, this is the same as:
  // ellK = ellf(PI/2.0,sqrt(k));
  ellK = rf(0.0,1.0-zmzp,1.0);
  
  // compute coordinate, proper, and Mino frequencies
  compute_freqs();
  
  // We will later evaluate many integrals of the form
  // int_0^x dx Xfunc(x).  It's easier to evaluate these
  // integrals using Chebyshev approximation of the integrand.
  prepare_Xfunc();

  // prepare coefficients for evaluating r(w_r) and theta(w_theta)
  //prepare_rtheta(); // TEST
  
  // prepare coefficients for evaluating non-oscilating parts of t and phi (Deltat and Deltaphi)
  prepare_Deltas();

  // compute r_isco and Omega_isco
  compute_isco();

  // compute average of cos^2 theta and of cot^2 theta
  compute_avgcos2();
  compute_avgcot2();
  
}
  
// overloaded constructor: this one lets you use E, Lz, and Q to specify the orbit
IEKG::IEKG(const Real spin, const Real E_in, const Real Lz_in, const Real Q_in, const Real tol) :  
	   a(spin), eps(tol), E(E_in), Lz(Lz_in), Q(Q_in)
{

  // compute cos(iota) as in [2]
  CosIota = Lz/sqrt(Lz*Lz + Q);

  // set prograde varriable
  if (Lz > 0) {
    prograde = 1;
  } else {
    prograde = 0;
  }

  // some functions in [1] needed to compute frequencies
  Real A = a*a*(1.0 - E*E);
  Real B = - Lz*Lz - Q - a*a*(1.0 - E*E);
  zPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*A);
  // Should be careful with the other root, which we'd like to define as
  // zMinusSquared = ( -B - sqrt(B*B - 4.0*A*Q) ) / (2.0*A); 
  // If AQ is small, ( -B - sqrt(B*B - 4.0*A*Q) ) is near zero (since B < 0 here)
  // see numerical recipies "Evaluation of Functions" for the following solution
  Real qq = -0.5*(B - sqrt(B*B - 4.0*A*Q) );
  zMinusSquared = Q/qq;
  zmzp = zMinusSquared / zPlusSquared;
  beta = a*sqrt(1.0 - E*E);
  Real betaOvera = sqrt(1.0 - E*E);
  aazPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*(1.0 - E*E));
  betaRootzPlusSquared = betaOvera*sqrt(aazPlusSquared);

  // define z_ as in [1] -- NOT AS IN [2,3]
  z_ = sqrt(zMinusSquared);

  // compute minimum theta
  theta_ = acos(z_);

  // compute rmin, rmax, e, and p 
  e_and_p();
  if (e < 0 || e > 1)
    Die("IEKG:: Eccentricity out of bounds.  Orbit has become unbound.");
  if (p < 0)
    Die("IEKG:: p < 0.  Something has gone wrong with the orbit...");

  // eccentricity, but cannot use an exactly zero eccentricity.
  if (e == 0)
    Die("IEKG:: Zero eccentricity?  Use the circular code, punky.");
	
  // some integrals in [1] needed to compute frequencies
  // Note, this is the same as:
  // ellK = ellf(PI/2.0,sqrt(k));
  ellK = rf(0.0,1.0-zmzp,1.0);
  
  // compute coordinate, proper, and Mino frequencies
  compute_freqs();
  
  // We will later evaluate many integrals of the form
  // int_0^x dx Xfunc(x).  It's easier to evaluate these
  // integrals using Chebyshev approximation of the integrand.
  prepare_Xfunc();

  // prepare coefficients for evaluating non-oscilating parts of t and phi (Deltat and Deltaphi)
  prepare_Deltas();

  // compute r_isco and Omega_isco
  compute_isco();

  // compute average of cos^2 theta and of cot^2 theta
  compute_avgcos2();
  compute_avgcot2();

}

// destructor
IEKG::~IEKG(){
  //free_Realvector(theta_k, 0, theta_kmax);
  //free_Realvector(r_n, 0, r_nmax);
  free_Realvector(Delta_t_k, 1, Delta_t_kmax);
  free_Realvector(Delta_t_n, 1, Delta_t_nmax);
  free_Realvector(Delta_phi_k, 1, Delta_phi_kmax);
  free_Realvector(Delta_phi_n, 1, Delta_phi_nmax);
  free_Realvector(cint,0,CN-1);
}

// compute <cos^2 theta>
void IEKG::compute_avgcos2(){
  if(theta_ == 0.5*M_PI) {
    avgcos2 = 0.0;
  } else {
    Real h_guess = M_PI / 100.0;
    avgcos2 = odeint(&IEKG::ignd_avgcos2, 0.0, M_PI, h_guess);
    avgcos2 /= M_PI;
  }
}
Real IEKG::ignd_avgcos2(Real chi){
  Real theta = Theta(chi);
  Real costheta = cos(theta);
  Real output = dwdchi(chi) * costheta * costheta;
  return output;
}

// compute <cot^2 theta>
void IEKG::compute_avgcot2(){
  if(theta_ == 0.5*M_PI) {
    avgcot2 = 0.0;
  } else {
    Real h_guess = M_PI / 100.0;
    avgcot2 = odeint(&IEKG::ignd_avgcot2, 0.0, M_PI, h_guess);
    avgcot2 /= M_PI;
  }
}
Real IEKG::ignd_avgcot2(Real chi){
  Real theta = Theta(chi);
  Real cottheta = 1.0/tan(theta);
  Real output = dwdchi(chi) * cottheta * cottheta;
  return output;
}

// compute r_isco and Omega_isco (see Finn & Thorne)
void IEKG::compute_isco(){
  Real Z1 = 1.0 + pow(1.0 - a*a,1.0/3.0) * ( pow(1.0 + a,1.0/3.0) + pow(1.0 - a,1.0/3.0) );
  Real Z2 = sqrt(3.0*a*a + Z1*Z1);
  if(prograde){
    r_isco = 3.0 + Z2 - sqrt( (3.0-Z1)*(3.0+Z1+2.0*Z2)  );
  }else{
    r_isco = 3.0 + Z2 + sqrt( (3.0-Z1)*(3.0+Z1+2.0*Z2)  );
  }
  Omega_isco = 1.0 / (pow(r_isco,1.5) + a);
}

// compute coordinate, proper, and Mino frequencies
void IEKG::compute_freqs(){
  // some integrals in [1] needed to compute frequencies
  // Note, this is the same as:
  // ellE = elle(PI/2.0,sqrt(k));
  // ellPi = ellpi(PI/2.0,-zMinusSquared,sqrt(k));
  ellK = rf(0.0,1.0-zmzp,1.0);
  Real ellE = ellK - zmzp*rd(0.0,1.0-zmzp,1.0)/3.0;
  Real ellPi = ellK + zMinusSquared*rj(0.0,1.0-zmzp,
                                       1.0,1.0-zMinusSquared)/3.0;

  // main batch of integrals from [1] 
  //
  // Romberg integration
  //Real W = qromb(&IEKG::Wfunc, 0.0, M_PI);
  //Real X = qromb(&IEKG::Xfunc, 0.0, M_PI);
  //Real Y = qromb(&IEKG::Yfunc, 0.0, M_PI);
  //Real Z = qromb(&IEKG::Zfunc, 0.0, M_PI);
  //
  // Adaptive stepsize integration
  Real h_guess = M_PI/1000.0;
  Real W = odeint(&IEKG::Wfunc, 0.0, M_PI, h_guess);
  Real X = odeint(&IEKG::Xfunc, 0.0, M_PI, h_guess);
  Real Y = odeint(&IEKG::Yfunc, 0.0, M_PI, h_guess);
  Real Z = odeint(&IEKG::Zfunc, 0.0, M_PI, h_guess);
  
  // another function from [1] needed to compute frequencies
  //Real Lambda = (Y + a*a*zPlusSquared*X)*ellK;
  //Lambda -= a*a*zPlusSquared*X*ellE;
  Real Lambda = (Y + aazPlusSquared*X)*ellK;
  Lambda -= aazPlusSquared*X*ellE;

  // action-angle frequencies
  OmegaR = M_PI*p*ellK / ((1.0 - e*e)*Lambda);
  //OmegaTheta = M_PI*beta*sqrt(zPlusSquared)*X / (2.0*Lambda);
  OmegaTheta = M_PI*betaRootzPlusSquared*X / (2.0*Lambda);
  OmegaPhi = ((Z - Lz*X)*ellK + Lz*X*ellPi  ) / Lambda;
  
  // Boyer-Lindquist frequencies
  //gamma = (W + a*a*zPlusSquared*E*X)*ellK;
  gamma = (W + aazPlusSquared*E*X)*ellK;
  //gamma -= a*a*zPlusSquared*E*X*ellE;
  gamma -= aazPlusSquared*E*X*ellE;
  gamma /= Lambda;
  OmegaR /= gamma;
  OmegaTheta /= gamma;
  OmegaPhi /= gamma;
  
  // Mino time frequencies as in [3]
  // NOTE: beta_here = sqrt(beta_[3]) and zPlusSquared_here = zPlus_[3].
  //UpsilonTheta = 0.5*M_PI*beta*sqrt(zPlusSquared) / ellK;
  UpsilonTheta = 0.5*M_PI*betaRootzPlusSquared / ellK;
  Gamma = UpsilonTheta / OmegaTheta;
  UpsilonR = Gamma*OmegaR;
  UpsilonPhi = Gamma*OmegaPhi;
}

// prepare for Chebyshev integration of Xfunc
void IEKG::prepare_Xfunc(){

  int MinCN = 40;
  int MaxCN = 1000;
  Real MinPsi = 0.0;
  Real MaxPsi = 2.0*M_PI;
  Real PsiCheck = 0.7*MaxPsi;

  // allocate some temporary memory for preperation
  Real *c; c=Realvector(0,MaxCN-1);
  Real *cintX; cintX=Realvector(0,MaxCN-1);

  // compute integrand using MinCN terms in the expansion
  CN = MinCN;
  chebft(MinPsi, MaxPsi, c, CN, &IEKG::Xfunc);
  chint(MinPsi, MaxPsi, c, cintX, CN);
  Real CheckInt = chebev(MinPsi, MaxPsi, cintX, CN, PsiCheck);
  Real FractionalChange = 2.0*eps;

  // establish accuracy in Chebyshev integration of Xfunc
  Real TempInt;
  while(FractionalChange > eps && CN < MaxCN-1){

    // increment CN
    CN++;

    // compute the new value of the integral
    chebft(MinPsi, MaxPsi, c, CN, &IEKG::Xfunc);
    chint(MinPsi, MaxPsi, c, cintX, CN);
    TempInt = chebev(MinPsi, MaxPsi, cintX, CN, PsiCheck);

    // compute the change and store new value
    FractionalChange = fabs(CheckInt - TempInt) / fabs(CheckInt);
    CheckInt = TempInt;

  }
  if(CN==MaxCN) Die("IEKG::  Too many terms in Chebyshev integration of Xfunc");

  // Fill the permanent array of Chebyshev integration coefs.
  cint = Realvector(0,CN-1);
  for(int j=0; j<CN; j++) cint[j] = cintX[j];

  // free the temporary memory from preping Chebyshev integration
  free_Realvector(c,0,MaxCN-1);
  free_Realvector(cintX,0,MaxCN-1);
}

// prepare coefficients for evaluating r(w_r) and theta(w_theta)
void IEKG::prepare_rtheta(){
  // make long array for temporary storage of expansion coeffecients
  // 100 or more terms needed likely means a problem.  
  int MaxCount = 200;  
  Real *theta_kX, *r_nX;
  theta_kX = Realvector(0,MaxCount);
  r_nX = Realvector(0,MaxCount);

  // Initialize the buffered truncation scheme.  The value of 
  // BufferLength is the number of consecutive terms which we require to 
  // make significantly small fractional contributions to the sums before
  // they are truncated.
  int BufferIndex = 0; 
  int BufferLength = 15;
  Real TruncatedSum = 0.0;
  Real FractionalChange = 2.0*eps;

  // If the orbit is exactly equatorial, theta_k vanishes for k != 0.
  if(theta_ == 0.5*M_PI){
      theta_kmax = 0;
      theta_kX[0] = theta_kcoef(theta_kmax);
  } else {
    // determine where to truncate theta(w_theta) angular expansion
    TruncatedSum = 0.0;
    FractionalChange = 2.0*eps;
    theta_kmax = -1;
    while(FractionalChange > eps){ 
      theta_kmax++;
      if(theta_kmax > MaxCount){
        //Die("Expansion sum in IEKG is too long: theta_kmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "theta_kmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if inclination is near zero."
             << "  If that's the case, ignore this warning.\n";
        theta_kmax = MaxCount;
        break;
      }
      theta_kX[theta_kmax] = theta_kcoef(theta_kmax);
      if(theta_kmax == 0){
        TruncatedSum += theta_kX[theta_kmax];
      } else {
        TruncatedSum += 2.0*theta_kX[theta_kmax];
      }
      FractionalChange = fabs(2.0*theta_kX[theta_kmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for 
      // truncation, we check that we have [BufferLength] terms in a row which 
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the equatorial limit, all terms in the sum vanish.  
      // We force this if after 100 terms, the sum is smaller than eps.
      if( theta_kmax==100 && fabs(TruncatedSum) < eps) {
        theta_kmax = 0;
        FractionalChange = 0.0;
      }
    }
  }

  // If the orbit is exactly circular, r_n vanish for n != 0.
  if(e < 1e-2){
      r_nmax = 0;
      r_nX[0] = r_ncoef(r_nmax);
  } else {
    // determine where to truncate r(w_r) expansion
    TruncatedSum = 0.0;
    FractionalChange = 1.0;
    r_nmax = -1;
    while(FractionalChange > eps){
      r_nmax++;
      if(r_nmax > MaxCount) {
        //Die("Expansion sum in IEKG is too long: Delta_t_nmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "r_nmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if eccentricity is near zero."
             << "  If that's the case, ignore this warning.\n";
        r_nmax = MaxCount;
        break;
      }
      r_nX[r_nmax] = r_ncoef(r_nmax);
      if(r_nmax == 0){
        TruncatedSum += r_nX[r_nmax];
      } else {
        TruncatedSum += 2.0*r_nX[r_nmax];
      }
      FractionalChange = fabs(2.0*r_nX[r_nmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for
      // truncation, we check that we have [BufferLength] terms in a row which
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the circular limit, all terms in the sum vanish.
      // We force this if after 100 terms, the sum is smaller than eps.
      if( r_nmax==100 && fabs(TruncatedSum) < eps) {
        r_nmax = 0;
        FractionalChange = 0.0;
      }
    } 
  }

  // Fill the permanent arrays of expansion coefficients
  theta_k = Realvector(0,theta_kmax);
  r_n = Realvector(0,r_nmax);
  for(int kk=0; kk<=theta_kmax; kk++){
    theta_k[kk] = theta_kX[kk];
    // TEST
    //cout.precision(16);
    //cout << theta_k[kk] << endl;    
  }
  for(int nn=0; nn<=r_nmax; nn++){
    r_n[nn] = r_nX[nn];
    // TEST
    //cout.precision(16);
    //cout << r_n[nn] << endl;
  }

  // Free the temporary long arrays of expansion coefficients
  free_Realvector(theta_kX, 0, MaxCount);
  free_Realvector(r_nX, 0, MaxCount);

}

// prepare coefficients for evaluating non-oscilating parts of t and phi (Deltat and Deltaphi)
void IEKG::prepare_Deltas(){
  // make long array for temporary storage of expansion coeffecients
  // 100 or more terms needed likely means a problem.  
  int MaxCount = 200;  
  Real *Delta_t_kX, *Delta_t_nX, *Delta_phi_kX, *Delta_phi_nX;
  Delta_t_kX = Realvector(1,MaxCount);
  Delta_t_nX = Realvector(1,MaxCount);
  Delta_phi_kX = Realvector(1,MaxCount);
  Delta_phi_nX = Realvector(1,MaxCount);

  // Initialize the buffered truncation scheme.  The value of 
  // BufferLength is the number of consecutive terms which we require to 
  // make significantly small fractional contributions to the sums before
  // they are truncated.
  int BufferIndex = 0; 
  int BufferLength = 15;
  Real TruncatedSum = 0.0;
  Real FractionalChange = 2.0*eps;

  // If the orbit is exactly equatorial, Delta_t^\theta vanishes.
  if(theta_ == 0.5*M_PI){
      Delta_t_kmax = 1;
      Delta_t_kX[1] = 0.0;
  } else {
    // determine where to truncate Delta_t angular expansion
    TruncatedSum = 0.0;
    FractionalChange = 2.0*eps;
    Delta_t_kmax = 0;
    while(FractionalChange > eps){ 
      Delta_t_kmax++;
      if(Delta_t_kmax > MaxCount){
        //Die("Expansion sum in IEKG is too long: Delta_t_kmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "Delta_t_kmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if inclination is near zero."
             << "  If that's the case, ignore this warning.\n";
        Delta_t_kmax = MaxCount;
        break;
      }
      Delta_t_kX[Delta_t_kmax] = Delta_t_kcoef(Delta_t_kmax);
      TruncatedSum += Delta_t_kX[Delta_t_kmax];
      FractionalChange = fabs(Delta_t_kX[Delta_t_kmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for 
      // truncation, we check that we have [BufferLength] terms in a row which 
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the equatorial limit, all terms in the sum vanish.  
      // We force this if after 100 terms, the sum is smaller than eps.
      if( Delta_t_kmax==100 && fabs(TruncatedSum) < eps) {
        Delta_t_kmax = 1;
        Delta_t_kX[1] = 0.0;
        FractionalChange = 0.0;
      }
    }
  } 

  // If the orbit is exactly equatorial, Delta_phi^\theta vanishes.
  if(theta_ == 0.5*M_PI){
      Delta_phi_kmax = 1;
      Delta_phi_kX[1] = 0.0;
  } else {
    // determine where to truncate Delta_phi angular expansion
    TruncatedSum = 0.0;
    FractionalChange = 1.0;
    Delta_phi_kmax = 0;
    while(FractionalChange > eps){
      Delta_phi_kmax++;
      if(Delta_phi_kmax > MaxCount) {
        //Die("Expansion sum in IEKG is too long:Delta_phi_kmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "Delta_phi_kmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if inclination is near zero."
             << "  If that's the case, ignore this warning.\n";
        Delta_phi_kmax = MaxCount;
        break;
      }
      Delta_phi_kX[Delta_phi_kmax] = Delta_phi_kcoef(Delta_phi_kmax);
      TruncatedSum += Delta_phi_kX[Delta_phi_kmax];
      FractionalChange = fabs(Delta_phi_kX[Delta_phi_kmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for
      // truncation, we check that we have [BufferLength] terms in a row which
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the equatorial limit, all terms in the sum vanish.
      // We force this if after 100 terms, the sum is smaller than eps.
      if( Delta_phi_kmax==100 && fabs(TruncatedSum) < eps) {
        Delta_phi_kmax = 1;
        Delta_phi_kX[1] = 0.0;
        FractionalChange = 0.0;
      }
    }
  }

  // If the orbit is exactly circular, Delta_t^r vanishes.
  if(e < 1e-2){
      Delta_t_nmax = 1;
      Delta_t_nX[1] = 0.0;
  } else {
    // determine where to truncate Delta_t radial expansion
    TruncatedSum = 0.0;
    FractionalChange = 1.0;
    Delta_t_nmax = 0;
    while(FractionalChange > eps){
      Delta_t_nmax++;
      if(Delta_t_nmax > MaxCount) {
        //Die("Expansion sum in IEKG is too long: Delta_t_nmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "Delta_t_nmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if eccentricity is near zero."
             << "  If that's the case, ignore this warning.\n";
        Delta_t_nmax = MaxCount;
        break;
      }
      Delta_t_nX[Delta_t_nmax] = Delta_t_ncoef(Delta_t_nmax);
      TruncatedSum += Delta_t_nX[Delta_t_nmax];
      FractionalChange = fabs(Delta_t_nX[Delta_t_nmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for
      // truncation, we check that we have [BufferLength] terms in a row which
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the circular limit, all terms in the sum vanish.
      // We force this if after 100 terms, the sum is smaller than eps.
      if( Delta_t_nmax==100 && fabs(TruncatedSum) < eps) {
        Delta_t_nmax = 1;
        Delta_t_nX[1] = 0.0;
        FractionalChange = 0.0;
      }
    }
  }

  // If the orbit is exactly circular, Delta_phi^r vanishes.
  if(e < 1e-2){
      Delta_phi_nmax = 1;
      Delta_phi_nX[1] = 0.0;
  } else {
    // determine where to truncate Delta_phi radial expansion
    TruncatedSum = 0.0;
    FractionalChange = 1.0;
    Delta_phi_nmax = 0;
    while (FractionalChange > eps){
      Delta_phi_nmax++;
      if(Delta_phi_nmax > MaxCount) {
        //Die("Expansion sum in IEKG is too long: Delta_phi_nmax > MaxCount");
        cout << "% -- WARNING -- Expansion sum in IEKG is too long: "
             << "Delta_phi_nmax > MaxCount = " << MaxCount << "\n"
             << "% This is to be expected if eccentricity is near zero."
             << "  If that's the case, ignore this warning.\n";
        Delta_phi_nmax = MaxCount;
        break;
      }
      Delta_phi_nX[Delta_phi_nmax] = Delta_phi_ncoef(Delta_phi_nmax);
      TruncatedSum += Delta_phi_nX[Delta_phi_nmax];
      FractionalChange = fabs(Delta_phi_nX[Delta_phi_nmax]/TruncatedSum);
      // buffer the truncation.  just to make sure the sum is really ready for
      // truncation, we check that we have [BufferLength] terms in a row which
      // give suitably small fractional contributions to the overall sum.
      if(FractionalChange > eps){
        BufferIndex=0;
      } else {
        BufferIndex++;
        FractionalChange = 2.0*eps;
        if(BufferIndex > BufferLength){
          FractionalChange = 0;
          BufferIndex = 0;
        }
      }
      // In the circular limit, all terms in the sum vanish.
      // We force this if after 100 terms, the sum is smaller than eps.
      if( Delta_phi_nmax==100 && fabs(TruncatedSum) < eps) {
        Delta_phi_nmax = 1;
        Delta_phi_nX[1] = 0.0;
        FractionalChange = 0.0;
      }
    }
  }

  // Fill the permanent arrays of expansion coefficients
  Delta_t_k = Realvector(1,Delta_t_kmax);
  Delta_t_n = Realvector(1,Delta_t_nmax);
  Delta_phi_k = Realvector(1,Delta_phi_kmax);
  Delta_phi_n = Realvector(1,Delta_phi_nmax);
  for(int kk=1; kk<=Delta_t_kmax; kk++){
    Delta_t_k[kk] = Delta_t_kX[kk];
  }
  for(int nn=1; nn<=Delta_t_nmax; nn++){
    Delta_t_n[nn] = Delta_t_nX[nn];
  }
  for(int kk=1; kk<=Delta_phi_kmax; kk++){
    Delta_phi_k[kk] = Delta_phi_kX[kk];
  }
  for(int nn=1; nn<=Delta_phi_nmax; nn++){
    Delta_phi_n[nn] = Delta_phi_nX[nn];
  }

  // Free the temporary long arrays of expansion coefficients
  free_Realvector(Delta_t_kX, 1, MaxCount);
  free_Realvector(Delta_t_nX, 1, MaxCount);
  free_Realvector(Delta_phi_kX, 1, MaxCount);
  free_Realvector(Delta_phi_nX, 1, MaxCount);
}

// m*Omega_phi + k*Omega_theta + n*Omega_r
Real IEKG::Omega(int mm, int kk, int nn)
{
  return mm*OmegaPhi + kk*OmegaTheta + nn*OmegaR;
}

// m*Upsilon_phi + k*Upsilon_theta + n*Upsilon_r
Real IEKG::Upsilon(int mm, int kk, int nn)
{
  return mm*UpsilonPhi + kk*UpsilonTheta + nn*UpsilonR;
}

// compute Delta_t(w_r,w_theta)
Real IEKG::Delta_t(Real w_r, Real w_theta)
{
  Real output = 0.0;
  for(int kk=1; kk<=Delta_t_kmax; kk++) output += Delta_t_k[kk]*sin(kk*w_theta);
  for(int nn=1; nn<=Delta_t_nmax; nn++) output += Delta_t_n[nn]*sin(nn*w_r);
  return output;
}
// overloaded version used to speed up calculations in IEKR 
Real IEKG::Delta_t(Real SgnRDot, Real SgnThetaDot, 
                   Real w_r, Real w_theta)
{
  Real r_term = 0.0;
  Real theta_term = 0.0;
  for(int kk=1; kk<=Delta_t_kmax; kk++) theta_term += Delta_t_k[kk]*sin(kk*w_theta);
  for(int nn=1; nn<=Delta_t_nmax; nn++) r_term += Delta_t_n[nn]*sin(nn*w_r);
  return SgnThetaDot*theta_term + SgnRDot*r_term;
}


// compute Delta_phi(chi,psi)
Real IEKG::Delta_phi(Real w_r, Real w_theta)
{
  Real output = 0.0;
  for(int kk=1; kk<=Delta_phi_kmax; kk++) output += Delta_phi_k[kk]*sin(kk*w_theta);
  for(int nn=1; nn<=Delta_phi_nmax; nn++) output += Delta_phi_n[nn]*sin(nn*w_r);
  return output;
}
// overloaded version used to speed up calculations in IEKR
Real IEKG::Delta_phi(Real SgnRDot, Real SgnThetaDot,
                     Real w_r, Real w_theta)
{
  Real theta_term = 0.0;
  Real r_term = 0.0;
  for(int kk=1; kk<=Delta_phi_kmax; kk++) theta_term += Delta_phi_k[kk]*sin(kk*w_theta);
  for(int nn=1; nn<=Delta_phi_nmax; nn++) r_term += Delta_phi_n[nn]*sin(nn*w_r);
  return SgnThetaDot*theta_term + SgnRDot*r_term;
}


// expansion coefficients
Real IEKG::theta_kcoef(int kk){
  k = kk; // set this in class so that integrand will see it
 
  Real output = 0.0;
  if(k==0){
    // theta_0 = <theta> = pi/2, by reflection symmetry: \dot \theta = f( \cos^2 \theta )
    output = 0.5*M_PI;
  } else {
    // use odeint
    Real h_guess = M_PI/(k*50.0);
    output = odeint(&IEKG::theta_k_integrand, 0.0, M_PI, h_guess);
    output /= M_PI;
  }
  return output;
}
Real IEKG::r_ncoef(int nn){
  n = nn; // set this in class so that integrand will see it

  // use odeint
  Real h_guess = (n == 0 ? M_PI / 50.0 : M_PI/(n*50.0) );
  Real output = odeint(&IEKG::r_n_integrand, 0.0, M_PI, h_guess);
  output /= M_PI;
  return output;
}
Real IEKG::Delta_t_kcoef(int kk)
{
  if(kk == 0){
    Die("IEKG::Delta_t_kcoef(k) called with argument k = 0.");
  }
  k = kk; // set this in class so that integrand will see it
  
  // use qromb
  //Real output = 2.0*qromb(&IEKG::Delta_t_k_integrand, 0.0, M_PI);
  //output /= k*M_PI*UpsilonTheta;
  //return output;

  // use odeint
  Real h_guess = 2.0*M_PI/(k*100.0);
  Real output = 2.0*odeint(&IEKG::Delta_t_k_integrand, 0.0, M_PI, h_guess);
  output /= k*M_PI*UpsilonTheta;
  return output;
}
Real IEKG::Delta_t_ncoef(int nn)
{

  if(nn == 0){
    Die("IEKG::Delta_t_ncoef(n) called with argument n = 0.");
  }
  n = nn; // set this in class so that integrand will see it

  // use qromb
  //Real output = 2.0*qromb(&IEKG::Delta_t_n_integrand, 0.0, M_PI);
  //output /= n*M_PI*UpsilonR;
  //return output;

  // use odeint
  Real h_guess = 2.0*M_PI/(n*100.0);
  Real output = 2.0*odeint(&IEKG::Delta_t_n_integrand, 0.0, M_PI,h_guess);
  output /= n*M_PI*UpsilonR;
  return output;

}
Real IEKG::Delta_phi_kcoef(int kk)
{
  if(kk == 0){
    Die("IEKG::Delta_phi_kcoef(k) called with argument k = 0.");
  }
  k = kk; // set this in class so that integrand will see it

  // use qromb
  //Real output = 2.0*qromb(&IEKG::Delta_phi_k_integrand, 0.0, M_PI);
  //output /= k*M_PI*UpsilonTheta;
  //return output;

  // use odeint
  Real h_guess = (k == 0 ? M_PI/50.0 : 2.0*M_PI/(k*100.0) );
  Real output = 2.0*odeint(&IEKG::Delta_phi_k_integrand, 0.0, M_PI, h_guess);
  output /= k*M_PI*UpsilonTheta;
  return output;
}
Real IEKG::Delta_phi_ncoef(int nn)
{
  if(nn == 0){
    Die("IEKG::Delta_phi_ncoef(n) called with argument n = 0.");
  }
  n = nn; // set this in class so that integrand will see it

  // use qromb
  //Real output = 2.0*qromb(&IEKG::Delta_phi_n_integrand, 0.0, M_PI);
  //output /= n*M_PI*UpsilonR;
  //return output;

  // use odeint
  Real h_guess = 2.0*M_PI/(n*100.0);
  Real output = 2.0*odeint(&IEKG::Delta_phi_n_integrand, 0.0, M_PI,h_guess);
  output /= n*M_PI*UpsilonR;
  return output;

}

//////////////////////////////////////////////////////////////////////////

// functions for integrals in [1]
Real IEKG::F(Real x)
{
  Real output = a*a*E/(p*p) - 2.0*a*(Lz-a*E)*(1.0 + e*cos(x))/(p*p*p);
  output *= (1.0 + e*cos(x))*(1.0 + e*cos(x));
  output += E;
  return output;
}
Real IEKG::G(Real x)
{
  return Lz - 2.0*(Lz - a*E)*(1.0+e*cos(x))/p;
}
Real IEKG::H(Real x)
{
  return 1.0 - 2.0*(1.0 + e*cos(x))/p +
    a*a*(1.0 + e*cos(x))*(1.0 + e*cos(x))/(p*p);
}
Real IEKG::J(Real x)
{
  Real output = (1.0-E*E)*(3.0+e*e)/(1.0-e*e) - 4.0/p +
    (a*a*(1.0-E*E) + Lz*Lz + Q)*(1.0-e*e)/(p*p);
  output *=(1.0+e*cos(x))*(1.0+e*cos(x));
  output += 2.0*(1.0 - E*E - (1.0-e*e)/p)*(1.0+e*cos(x));
  output += (1.0-E*E)*(1.0-e*e);
  return output;
}
Real IEKG::Wfunc(Real x)
{
  return p*p*F(x)/((1.0 + e*cos(x))*(1.0 + e*cos(x))*H(x)*
		   sqrt(J(x)));
}
Real IEKG::Xfunc(Real x)
{
  return 1.0/sqrt(J(x));
}
Real IEKG::Yfunc(Real x)
{
  return p*p/((1.0 + e*cos(x))*(1.0 + e*cos(x))*sqrt(J(x)));
}
Real IEKG::Zfunc(Real x)
{
  return G(x) / (H(x)*sqrt(J(x)));
}

// deritvative of w_theta with respect to chi
Real IEKG::dwdchi(Real chi)
{
  //Real output = 1.0/sqrt(1.0 - zMinusSquared*cos(chi)*cos(chi)/zPlusSquared);
  Real output = 1.0/sqrt(1.0 - zmzp*cos(chi)*cos(chi));
  output *= 0.5*M_PI / ellK;
  return output;
}

// deritvative of w_r with respect to psi
Real IEKG::dwdpsi(Real psi)
{
  return UpsilonR*(1.0 - e*e) / ( p*sqrt(J(psi)) );
}

// d t/d lambda(chi,psi) = TofChi(chi) + TofPsi(psi)
Real IEKG::TofChi(Real chi)
{
  Real theta = Theta(chi);

  return -E*a*a*sin(theta)*sin(theta);
}
Real IEKG::TofPsi(Real psi)
{
  Real r = RofPsi(psi);
  Real Delta = r*r - 2.0*r + a*a;
  return E*(r*r + a*a)*(r*r + a*a)/Delta
         + a*Lz*(1.0 - (r*r + a*a)/Delta );
}

// d phi/d lambda(chi,psi) = PhiofChi(chi) + PhiofPsi(psi)
Real IEKG::PhiofChi(Real chi)
{
  Real theta = Theta(chi);
  return Lz /( sin(theta)*sin(theta) );
}
Real IEKG::PhiofPsi(Real psi)
{
  Real r = RofPsi(psi);
  Real Delta = r*r - 2.0*r + a*a;
  return a*E*( (r*r + a*a)/Delta  - 1.0) - a*a*Lz/Delta;
}

// w_theta as a function of 0 < chi < 2pi
Real IEKG::wTheta(Real chi)
{
  Real lambda;
  if(chi >= 0 && chi <= 0.5*M_PI){
    lambda = lambda0(chi);
  } else if(chi > 0.5*M_PI && chi <= M_PI){
    //lambda = 2.0*ellK/(beta*sqrt(zPlusSquared)) - lambda0(M_PI - chi);
    lambda = 2.0*ellK/betaRootzPlusSquared - lambda0(M_PI - chi);
  } else if(chi > M_PI && chi <= 1.5*M_PI){
    //lambda = 2.0*ellK/(beta*sqrt(zPlusSquared)) + lambda0(chi - M_PI);
    lambda = 2.0*ellK/betaRootzPlusSquared + lambda0(chi - M_PI);
  } else if(chi > 1.5*M_PI && chi <= 2.0*M_PI){
    //lambda = 4.0*ellK/(beta*sqrt(zPlusSquared)) - lambda0(2.0*M_PI - chi);
    lambda = 4.0*ellK/betaRootzPlusSquared - lambda0(2.0*M_PI - chi);
  } else {
    Die("chi out of range in wTheta(chi).  we require 0 < chi < 2pi.");
  }
  return UpsilonTheta*lambda;
}

// w_r as a function of psi
Real IEKG::wR(Real psi)
{

  // use 10-point Gaussian quatrature
  //Real lambda = (1.0 - e*e) * qgaus(&IEKG::Xfunc, 0.0, psi) / p;
  //return UpsilonR*lambda;

  // use qromb 
  //Real lambda = (1.0 - e*e) * qromb(&IEKG::Xfunc, 0.0, psi) / p;
  //return UpsilonR*lambda;

  // use odeint 
  //Real h_guess = psi/10.0;
  //Real lambda = (1.0 - e*e) * odeint(&IEKG::Xfunc, 0.0, psi, h_guess) / p;
  //return UpsilonR*lambda;

  // use Chebyshev integration
  Real lambda = (1.0 - e*e) * chebev(0.0, 2.0*M_PI, cint, CN, psi) / p;
  return UpsilonR*lambda;
}

// theta as a function of chi
// NOTE: zPlusSquared_here = zPlus_[3]
Real IEKG::Theta(Real chi)
{
  return acos( sqrt(zMinusSquared)*cos(chi) );
}

// r as a function of psi
Real IEKG::RofPsi(Real psi)
{
  return p / ( 1.0 + e*cos(psi) );
}

// lambda(chi) for 0 < chi < pi/2
// NOTE: beta_here = sqrt(beta_[3]) and zPlusSquared_here = zPlus_[3]
// and zMinusSquared = zMinus_[3]
Real IEKG::lambda0(Real chi)
{
  Real output = ellK;
  //output -= ellf(0.5*M_PI - chi, sqrt(zMinusSquared/zPlusSquared));
  output -= ellf(0.5*M_PI - chi, sqrt(zmzp));
  //output /= beta*sqrt(zPlusSquared); 
  output /= betaRootzPlusSquared;
  return output;
}

// integrands for expansion coefficients
Real IEKG::theta_k_integrand(Real chi){
  Real w_theta = wTheta(chi);
  return dwdchi(chi) * Theta(chi) * cos(k*w_theta);
}
Real IEKG::r_n_integrand(Real psi){
  Real w_r = wR(psi);
  return dwdpsi(psi) * RofPsi(psi) * cos(n*w_r);
}
Real IEKG::Delta_t_k_integrand(Real chi)
{
  Real w_theta = wTheta(chi);
  return dwdchi(chi) * TofChi(chi) * cos(k*w_theta);
}
Real IEKG::Delta_t_n_integrand(Real psi)
{
  Real w_r = wR(psi);
  return dwdpsi(psi) * TofPsi(psi) * cos(n*w_r);
}
Real IEKG::Delta_phi_k_integrand(Real chi)
{
  Real w_theta = wTheta(chi);
  return dwdchi(chi) * PhiofChi(chi) * cos(k*w_theta);
}
Real IEKG::Delta_phi_n_integrand(Real psi)
{
  Real w_r = wR(psi); 
  return dwdpsi(psi) * PhiofPsi(psi) * cos(n*w_r);
}

//
// Adaptive step size integration routines: odeint, rkqs, rkck
//
// from Numerical Recipes
//
Real IEKG::odeint(Real (IEKG::*derivs)(Real), Real x1, Real x2, Real h1)
{
        const int MAXSTP=1000000;
        const Real TINY=1.0e-30;
        Real hmin = 0;
        int i,nstp;
        Real x,hnext,hdid,h;

        Real yscal; 
        Real y; 
        Real dydx;
        x=x1;
        h=SIGN(h1,x2-x1);
        y = 0.0;
        for (nstp=0;nstp<MAXSTP;nstp++) {
                dydx=(this->*derivs)(x);
                yscal=fabs(y)+fabs(dydx*h)+TINY;
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs);
                if ((x-x2)*(x2-x1) >= 0.0) {
                        return y;
                }
                if (fabs(hnext) <= hmin) Die("Step size too small in IEKG::odeint");
                h=hnext;
        }
        Die("Too many steps in routine IEKG::odeint");
        return 0.0; // never make it here. this just quiets a compile warning
}
void IEKG::rkqs(Real &y, Real dydx, Real &x, Real htry,
        Real eps, Real yscal, Real &hdid, Real &hnext,
        Real (IEKG::*derivs)(Real) )
{
        const Real SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
        int i;
        Real errmax,h,htemp,xnew;

        h=htry;
        Real yerr; 
        Real ytemp;
        for (;;) {
                rkck(y,dydx,x,h,ytemp,yerr,derivs);
                errmax=0.0;
                errmax=Max(errmax,fabs(yerr/yscal));
                errmax /= eps;
                if (errmax <= 1.0) break;
                htemp=SAFETY*h*pow(errmax,PSHRNK);
                h=(h >= 0.0 ? Max(htemp,0.1*h) : Min(htemp,0.1*h));
                xnew=x+h;
                if (xnew == x) Die("stepsize underflow in IEKG::rkqs");
        }
        if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
        else hnext=5.0*h;
        x += (hdid=h);
        y=ytemp;
}
void IEKG::rkck(Real &y, Real &dydx, Real x,
        Real h, Real &yout, Real &yerr,
        Real (IEKG::*derivs)(Real) )
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

        Real ak2;
        Real ak3;
        Real ak4;
        Real ak5;
        Real ak6;
        ak2=(this->*derivs)(x+a2*h);
        ak3=(this->*derivs)(x+a3*h);
        ak4=(this->*derivs)(x+a4*h);
        ak5=(this->*derivs)(x+a5*h);
        ak6=(this->*derivs)(x+a6*h);
        yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6);
        yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);
}

//
// Romberg integration from Numerical Recipes (uses trpzd and polint)
//
// NOTE: this has been modified to use a trpzd function without
//       static varriables.
//
Real IEKG::qromb(Real (IEKG::*func)(Real), Real start, Real finish)
{
        int JMAX=100, JMAXP=JMAX+1, K=6;
        Real ss,dss;
        Real *s, *h, *s_t, *h_t;
        s = Realvector(0,JMAX-1);
        h = Realvector(0,JMAXP-1);
        s_t = Realvector(0,K-1);
        h_t = Realvector(0,K-1);
        int i,j;

        h[0]=1.0;
        for (j=1;j<=JMAX;j++) {
                if(j==1){
                  trpzd(func,start,finish,j,0.0,&s[0]);
                } else {
                  trpzd(func,start,finish,j,s[j-2],&s[j-1]);
                }
                if (j >= K) {
                        for (i=0;i<K;i++) {
                                h_t[i]=h[j-K+i];
                                s_t[i]=s[j-K+i];
                        }
                        polint(h_t,s_t,0.0,ss,dss,K);
                        //if (fabs(dss) <= eps*fabs(ss)) return ss;
                        //modified from above to account for 0 result
                        if (fabs(dss) <= eps*fabs(ss) || fabs(ss) < eps){
                          free_Realvector(s,0,JMAX-1);
                          free_Realvector(h,0,JMAXP-1);
                          free_Realvector(s_t,0,K-1);
                          free_Realvector(h_t,0,K-1);
                          return ss;
                        }
                }
                h[j]=0.25*h[j-1];
        }
        Die("Too many steps in routine qromb");
        return 0.0;
}
void IEKG::trpzd(Real (IEKG::*func)(Real), Real start, Real finish, int nn, Real last, Real *next)
{
        Real x,tnm,sum,del;
        int it,j;
                                                                                                                        
        if (nn == 1) {
                *next=0.5*(finish-start)*(   (this->*func)(start) + (this->*func)(finish)  );
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
void IEKG::polint(Real *xa, Real *ya, Real x, Real &y, Real &dy, int nn)
{
        int i,mm,ns=0;
        Real den,dif,dift,ho,hp,w;
        Real *c, *d;
        c = Realvector(0,nn-1);
        d = Realvector(0,nn-1);
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
                        if ((den=ho-hp) == 0.0) Die("Error in routine polint");
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                y += (dy=(2*(ns+1) < (nn-mm) ? c[ns+1] : d[ns--]));
        }
        free_Realvector(c,0,nn-1);
        free_Realvector(d,0,nn-1);
}
//
// Gaussian quadrature routine from Numerical Recipes
// (for the wild at heart)
//
Real IEKG::qgaus(Real (IEKG::*func)(Real), Real start, Real finish)
{
  static const Real x[]={0.1488743389816312, 0.4333953941292472,
                         0.6794095682990244, 0.8650633666889845,
                         0.9739065285171717};
  static const Real w[]={0.2955242247147529, 0.2692667193099963,
                         0.2190863625159821, 0.1494513491505806,
                         0.0666713443086881};
  int j;
  Real xr, xm, dx, s;

  xm = 0.5*(finish + start);
  xr = 0.5*(finish - start);
  s=0;
  for (j = 0; j < 5; j++) {
    dx = xr*x[j];
    s += w[j]*(  (this->*func)(xm + dx) + (this->*func)(xm - dx)  );
  }
  return s *= xr;
}

//
// Chebyshev Integration routines
//
void IEKG::chebft(Real aa, Real b, Real *c, int nn, Real (IEKG::*func)(Real))
{
        int kk,j;
        Real fac,bpa,bma,y,sum;

        Real *f; f=Realvector(0,nn);
        bma=0.5*(b-aa);
        bpa=0.5*(b+aa);
        for (kk=0;kk<nn;kk++) {
                y=cos(M_PI*(kk+0.5)/nn);
                f[kk]=(this->*func)(y*bma+bpa);
        }
        fac=2.0/nn;
        for (j=0;j<nn;j++) {
                sum=0.0;
                for (kk=0;kk<nn;kk++)
                        sum += f[kk]*cos(M_PI*j*(kk+0.5)/nn);
                c[j]=fac*sum;
        }
        free_Realvector(f,0,nn);
}
Real IEKG::chebev(Real aa, Real b, Real *c, int m, Real x)
{
        Real d=0.0,dd=0.0,sv,y,y2;
        int j;

        if ((x-aa)*(x-b) > 0.0)
                Die("IEKG::chebev: x not in range in routine");
        y2=2.0*(y=(2.0*x-aa-b)/(b-aa));
        for (j=m-1;j>0;j--) {
                sv=d;
                d=y2*d-dd+c[j];
                dd=sv;
        }
        return y*d-dd+0.5*c[0];
}
void IEKG::chint(Real aa, Real b, Real *c, Real *cint, int nn)
{
        int j;
        Real sum=0.0,fac=1.0,con;

        con=0.25*(b-aa);
        for (j=1;j<nn-1;j++) {
                cint[j]=con*(c[j-1]-c[j+1])/j;
                sum += fac*cint[j];
                fac = -fac;
        }
        cint[nn-1]=con*c[nn-2]/(nn-1);
        sum += fac*cint[nn-1];
        cint[0]=2.0*sum;
}

// determine e and p from E Lz and Q
// converted from Scott's "Rroots" in an early version of IEKG
// most of the comments and code here are Scott's.
void IEKG::e_and_p()
{
  int N_R = 4;
  Complex *R_cofs, *R_roots;
  Real ome2 = (1. - E)*(1. + E);
  //
  // Store coefficients of polynomial for R in array, find roots.
  //
  R_cofs = Complexvector(0, N_R);
  R_cofs[0] = -a*a*Q/ome2;
  R_cofs[1] = 2.*(Q + (a*E - Lz)*(a*E - Lz))/ome2;
  R_cofs[2] = -(a*a*ome2 + Lz*Lz + Q)/ome2;
  R_cofs[3] = 2./ome2;
  R_cofs[4] = -1.;
  R_roots = Complexvector(1, N_R);
  zroots(R_cofs, N_R, R_roots, 1, eps);
  free_Complexvector(R_cofs, 0, N_R);
  //
  // Remap the roots to what I've got in my notes.  Note that with
  // the current scheme for inputting parameters we NEVER get an
  // imaginary part to the roots (well, maybe a small one due to
  // roundoff error).
  //
  Real arr1 = R_roots[4].real();
  Real arr2 = R_roots[3].real();
  Real arr3 = R_roots[2].real();
  Real arr4 = R_roots[1].real();

  // store minimum and maximum radius
  rmax = arr1;
  rmin = arr2;
  
  // Store eccentricity and semi-latus rectum.
  // There are lots of different ways to do this. 
  // I'm not sure which evaluates best.
  //
  p = 2.0*rmax*rmin / (rmax + rmin);
  e = (rmax - rmin) / (rmax + rmin);
  //e = (1.0 - rmin/rmax) / (1.0 + rmin/rmax);
  //e = 1.0 - (p/rmax);
  //e = (p/rmin) - 1.0;
  //e = 0.5*p*( (1.0/rmin) - (1.0/rmax) );
  //p = rmin*(1+e);

  free_Complexvector(R_roots, 1, N_R);
}
void IEKG::zroots(Complex a[], int m, Complex roots[], int polish, Real EPS)
{
  int MAXM = 100;
  int i, its, j, jj;
  Complex x, b, c, ad[MAXM];

  for (j = 0; j <= m; j++)
    ad[j] = a[j];
  for (j = m; j >= 1; j--) {
    x = Complex(0.0, 0.0);
    laguer(ad, j, &x, &its);
    if (fabs(x.imag()) <= 2.0*EPS*fabs(x.real()))
      x = Complex(x.real(), 0.0);
    roots[j] = x;
    b = ad[j];
    for (jj = j - 1; jj >= 0; jj--) {
      c = ad[jj];
      ad[jj] = b;
      b = c + x*b;
    }
  }
  if (polish)
    for (j = 1; j <= m; j++)
      laguer(a, m, &roots[j], &its);
  for (j = 2; j <= m; j++) {
    x = roots[j];
    for (i = j - 1; i >= 1; i--) {
      if (roots[i].real() <= x.real()) break;
      roots[i + 1] = roots[i];
    }
    roots[i + 1] = x;
  }
}
void IEKG::laguer(Complex a[], int m, Complex *x, int *its)
{
  Real EPSS=1.0e-15;
  int MR=15;
  int MT=20;
  int MAXIT=(MT*MR);
  int iter, j;
  Real abx, abp, abm, err;
  Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
  static Real frac[16] =
  {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

  for (iter = 1; iter <= MAXIT; iter++) {
    *its = iter;
    b = a[m];
    err = abs(b);
    d = f = Complex(0.0, 0.0);
    abx = abs(*x);
    for (j = m - 1; j >= 0; j--) {
      f = (*x)*f + d;
      d = (*x)*d + b;
      b = (*x)*b + a[j];
      err = abs(b) + abx*err;
    }
    err *= EPSS;
    if (abs(b) <= err) return;
    g = d/b;
    g2 = g*g;
    h = g2 - (2.0 * f/b);
    sq = sqrt(((Real)(m - 1))*(((Real)m)*h - g2));
    gp = g + sq;
    gm = g - sq;
    abp = abs(gp);
    abm = abs(gm);
    if (abp < abm) gp = gm;
    dx = (FMAX(abp, abm) > 0.0 ?
          Complex((float) m, 0.0)/gp :
          exp(log(1 + abx)) * Complex(cos((float)iter), sin((float)iter)));
    x1 = (*x) - dx;
    if (x->real() == x1.real() && x->imag() == x1.imag())
      return;
    if (iter % MT)
      *x = x1;
    else
      *x = *x - frac[iter/MT]*dx;
  }
  Die("too many iterations in IEKG::laguer");
  return;
}
