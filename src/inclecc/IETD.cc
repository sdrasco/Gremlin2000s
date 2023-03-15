//
// Steve Drasco, 5 Jan 2006
//

// standard includes
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>
#include <math.h>

// Gremlin includes
#include "SNT.h"
#include "SWSH.h"
#include "Globals.h"
#include "IEKG.h"
#include "IEKR.h"
#include "IETD.h"

// constructor
IETD::IETD(IEKG *iekg_in, Real eps_in) : iekg(iekg_in), eps(eps_in)
{

  // define global varriables
  MAXSTEP = 200;
  USEnSlope = 1;
  USEkSlope = 1;
  K_BUFF = 3;
  N_BUFF = 5;
  TEPS = 0.001;
  EPS_L = eps;
  EPS_K = eps;
  EPS_N = eps;
  MaxEIlmkn = 0.0;
  MaxEIlmk  = 0.0;
  MaxEIlm   = 0.0;
  MaxEIl    = 0.0;
  MaxEHlmkn = 0.0;
  MaxEHlmk  = 0.0;
  MaxEHlm   = 0.0;
  MaxEHl    = 0.0;
  MaxLIlmkn = 0.0;
  MaxLIlmk  = 0.0;
  MaxLIlm   = 0.0;
  MaxLIl    = 0.0;
  MaxLHlmkn = 0.0;
  MaxLHlmk  = 0.0;
  MaxLHlm   = 0.0;
  MaxLHl    = 0.0;
  TempMaxEIlmkn = 0.0;
  kMax = 0;
  nMax = 0;

  dEH = 0.0;
  dEI = 0.0;
  dLH = 0.0;
  dLI = 0.0;
  dQ = 0.0;

  Xl = Realvector(0,3);
  Xlm = Realvector(0,3);
  Xlmk = Realvector(0,3);
  Xlmkn = Realvector(0,3);

  // set initial nMax using Peters-Mathews formula // TEST
  //if(iekg->e > 0.1) {
  //  nMax = int( exp( 0.5 - 1.5*log(1 - iekg->e)) );
  //}

  // open output file and print header information (as matlab comment)
  char fname[50];
  sprintf(fname, "IETD_output_%f_%f_%f_%f", iekg->a,iekg->e,iekg->p,(iekg->theta_)*180/M_PI);
  outfile.open(fname,ios::out);
  outfile.precision(16);
  outfile << "%\n"
       << "%a = " << iekg->a << "\n"
       << "%e = " << iekg->e << "\n"
       << "%p = " << iekg->p << "\n"
       << "%MinTheta = " << (iekg->theta_)*180/M_PI << " (degrees)\n"
       << "%IEKG.CosIota = " << iekg->CosIota << endl
       << "%Iota = " << acos(iekg->CosIota) * 180.0/M_PI <<  " (degrees)\n"
       << "%E = " << iekg->E << endl
       << "%Lz = " << iekg->Lz << endl
       << "%Q = " << iekg->Q << endl
       << "%OmegaR = " << iekg->OmegaR << endl
       << "%OmegaTheta = " << iekg->OmegaTheta << endl
       << "%OmegaPhi = " << iekg->OmegaPhi << endl
       << "%gamma = " << iekg->Gamma << endl
       << "%UpsilonR = " << iekg->UpsilonR << endl
       << "%UpsilonTheta = " << iekg->UpsilonTheta << endl
       << "%UpsilonPhi = " << iekg->UpsilonPhi << endl
       << "%eps = " << eps << "\n"
       << "%eps_l = " << EPS_L << "\n"
       << "%eps_k = " << EPS_K << "\n"
       << "%eps_n = " << EPS_N << "\n"
       << "%Teps = " << TEPS << "\n"
       << "%Delta_t_kmax = " << iekg->Delta_t_kmax << "\n"
       << "%Delta_t_nmax = " << iekg->Delta_t_nmax << "\n"
       << "%Delta_phi_kmax = " << iekg->Delta_phi_kmax << "\n"
       << "%Delta_phi_nmax = " << iekg->Delta_phi_nmax << "\n"
       << "%MAXSTEP = " << MAXSTEP << "\n"
       << "%K_BUFF = " << K_BUFF << "\n"
       << "%N_BUFF = " << N_BUFF << "\n"
       << "%USEnSlope = " << USEnSlope << "\n"
       << "%USEkSlope = " << USEkSlope << "\n"
       << "%\n"
       << "%below: [l m k n omega Z^H.real() ZH.imag() Z^Inf.real() ZInf.imag() ... \n"
//       << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH QInf QH]" << "\n" // TEST
       << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH]" << "\n"
       << "%\n";
  outfile.flush();

  // outermost loop over l-index is done here
  // data is recorded internaly
  Real epsl = EPS_L;
  int l = 2, ContinueSumH = 1, ContinueSumI = 1;
  while(ContinueSumH || ContinueSumI){

    // reset temporary peak value of (dE/dt)^{\infty}
    // for tracking peaks in n-index, since the peak's location 
    // looks linear in l-index.
    TempMaxEIlmkn = 0.0;

    // get next set of terms in sum
    GetXl(l, ContinueSumH, ContinueSumI);

    // check if we have any new maxima
    if(fabs(Xl[0]) > MaxEIl) MaxEIl = fabs(Xl[0]);
    if(fabs(Xl[1]) > MaxEHl) MaxEHl = fabs(Xl[1]);
    if(fabs(Xl[2]) > MaxLIl) MaxLIl = fabs(Xl[2]);
    if(fabs(Xl[3]) > MaxLHl) MaxLHl = fabs(Xl[3]);

    // increment l-index
    l++;

    // check truncation criteria
    ContinueSumH = (fabs(Xl[1]) > epsl*MaxEHl ||
                    fabs(Xl[3]) > epsl*MaxLHl   );
    ContinueSumI = (fabs(Xl[0]) > epsl*MaxEIl ||
                    fabs(Xl[2]) > epsl*MaxLIl   );

    // we only take so many steps
    if(l-2 > MAXSTEP) Die("MAXSTEP exceeded!");

  }

}

// destructor
IETD::~IETD()
{
  free_Realvector(Xl, 0, 3);
  free_Realvector(Xlm, 0, 3);
  free_Realvector(Xlmk, 0, 3);
  free_Realvector(Xlmkn, 0, 3);
  if (outfile.is_open()) {
    outfile.close();
  }
}

void IETD::GetXl(int l, int ContinueSumH, int ContinueSumI){

  // initialize stuff
  Xl[0]=Xl[1]=Xl[2]=Xl[3]=0.0;
  Real epsm = 100.0*eps;

  //
  // go through all possible values for m
  //
  for(int m=-l; m<=l; m++){
    // get new terms
    GetXlm(l,m,ContinueSumH,ContinueSumI);

    // increment output terms
    Xl[0] += Xlm[0];
    Xl[1] += Xlm[1];
    Xl[2] += Xlm[2];
    Xl[3] += Xlm[3];
  }

}

void IETD::GetXlm(int l, int m, int DoSumH, int DoSumI){

  // initialize stuff
  Real InitialX0, InitialX1, InitialX2, InitialX3;
  int BufferIndexH=0,BufferIndexI=0,BufferLength=K_BUFF;
  Xlm[0]=Xlm[1]=Xlm[2]=Xlm[3]=0.0;
  int kMaxOld = kMax;
  Real epsk = EPS_K;

  // forward loop (positive)
  int ContinueSumH=DoSumH, ContinueSumI=DoSumI;
  int k=kMax; // start at location of last maximum
  while(ContinueSumH || ContinueSumI){

    // get new terms
    GetXlmk(l,m,k,ContinueSumH,ContinueSumI);

    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmk[0]) > MaxEIlmk){
      MaxEIlmk = fabs(Xlmk[0]);
      kMax = k;
    }
    if(fabs(Xlmk[1]) > MaxEHlmk) MaxEHlmk = fabs(Xlmk[1]);
    if(fabs(Xlmk[2]) > MaxLIlmk) MaxLIlmk = fabs(Xlmk[2]);
    if(fabs(Xlmk[3]) > MaxLHlmk) MaxLHlmk = fabs(Xlmk[3]);

    // increment k-index
    k++;

    // check truncation criteria
    ContinueSumH = ( fabs(Xlmk[1]) > epsk*MaxEHlmk ||
                     fabs(Xlmk[3]) > epsk*MaxLHlmk   );
    ContinueSumI = ( fabs(Xlmk[0]) > epsk*MaxEIlmk ||
                     fabs(Xlmk[2]) > epsk*MaxLIlmk   );

    // buffer the truncation (H)
    if(DoSumH && BufferIndexH <= BufferLength){
      if(ContinueSumH){
        BufferIndexH=0;
      } else {
        if(BufferIndexH==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1 = fabs(Xlmk[1]); // EH
          InitialX3 = ( m==0 ? 1.0 : fabs(Xlmk[3]) ); // LH
        }
        BufferIndexH++;
        ContinueSumH=DoSumH;
        if(BufferIndexH > BufferLength){
          ContinueSumH = 0;
          if((InitialX1 < fabs(Xlmk[1]) ||  // EH modes are increasing
              InitialX3 < fabs(Xlmk[3]) )   // LH modes are increasing
            && USEkSlope  ) {
            ContinueSumH = 1;
            BufferIndexH = 0;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(DoSumI && BufferIndexI <= BufferLength){
      if(ContinueSumI){
        BufferIndexI=0;
      } else {
        if(BufferIndexI==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0 = fabs(Xlmk[0]); // EI
          InitialX2 = ( m==0 ? 1.0 : fabs(Xlmk[2]) ); // LI
        }
        BufferIndexI++;
        ContinueSumI=DoSumI;
        if(BufferIndexI > BufferLength){
          ContinueSumI = 0;
          if((InitialX0 < fabs(Xlmk[0]) ||  // EI modes are increasing
              InitialX2 < fabs(Xlmk[2]) )   // LI modes are increasing
            && USEkSlope  ) {
            ContinueSumI = 1;
            BufferIndexI = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(k-kMax) > MAXSTEP) Die("MAXSTEP exceeded!");

    // increment output terms
    Xlm[0] += Xlmk[0];
    Xlm[1] += Xlmk[1];
    Xlm[2] += Xlmk[2];
    Xlm[3] += Xlmk[3];

  }

  // backward loop (negative)
  k=kMaxOld-1;  // start at location of last maximum
  ContinueSumH=DoSumH, ContinueSumI=DoSumI; // make sure the loop starts
  BufferIndexH = 0, BufferIndexI = 0;
  while(ContinueSumH || ContinueSumI){

    // get new terms
    GetXlmk(l,m,k,ContinueSumH,ContinueSumI);

    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmk[0]) > MaxEIlmk){
      MaxEIlmk = fabs(Xlmk[0]);
      kMax = k;
    }
    if(fabs(Xlmk[1]) > MaxEHlmk) MaxEHlmk = fabs(Xlmk[1]);
    if(fabs(Xlmk[2]) > MaxLIlmk) MaxLIlmk = fabs(Xlmk[2]);
    if(fabs(Xlmk[3]) > MaxLHlmk) MaxLHlmk = fabs(Xlmk[3]);

    // increment k-index
    k--;

    // check truncation criteria
    ContinueSumH = ( fabs(Xlmk[1]) > epsk*MaxEHlmk ||
                     fabs(Xlmk[3]) > epsk*MaxLHlmk   );
    ContinueSumI = ( fabs(Xlmk[0]) > epsk*MaxEIlmk ||
                    fabs(Xlmk[2]) > epsk*MaxLIlmk   );

    // buffer the truncation (H)
    if(DoSumH && BufferIndexH <= BufferLength){
      if(ContinueSumH){
        BufferIndexH=0;
      } else {
        if(BufferIndexH==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1 = fabs(Xlmk[1]); // EH
          InitialX3 = ( m==0 ? 1.0 : fabs(Xlmk[3]) ); // LH
        }
        BufferIndexH++;
        ContinueSumH=1;
        if(BufferIndexH > BufferLength){
          ContinueSumH = 0;
          if((InitialX1 < fabs(Xlmk[1]) ||  // EH modes are increasing
              InitialX3 < fabs(Xlmk[3]) )   // LH modes are increasing
            && USEkSlope  ) {
            ContinueSumH = 1;
            BufferIndexH = 0;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(DoSumI && BufferIndexI <= BufferLength){
      if(ContinueSumI){
        BufferIndexI=0;
      } else {
        if(BufferIndexI==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0 = fabs(Xlmk[0]); // EI
          InitialX2 = ( m==0 ? 1.0 : fabs(Xlmk[2]) ); // LI
        }
        BufferIndexI++;
        ContinueSumI=1;
        if(BufferIndexI > BufferLength){
          ContinueSumI = 0;
          if((InitialX0 < fabs(Xlmk[0]) ||  // EI modes are increasing
              InitialX2 < fabs(Xlmk[2]) )   // LI modes are increasing
            && USEkSlope  ) {
            ContinueSumI = 1;
            BufferIndexI = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(k-kMax) > MAXSTEP) Die("MAXSTEP exceeded!");

    // increment output terms
    Xlm[0] += Xlmk[0];
    Xlm[1] += Xlmk[1];
    Xlm[2] += Xlmk[2];
    Xlm[3] += Xlmk[3];

  }

}

void IETD::GetXlmk(int l, int m, int k, int DoSumH, int DoSumI){

  // if this is an equatorial orbit, only k=0 terms are nonzero. 
  if(iekg->theta_ == 0.5*M_PI && k != 0) {
    Xlmk[0]=Xlmk[1]=Xlmk[2]=Xlmk[3]=0.0;
    return;
  }

  // initialize stuff
  Real InitialX0, InitialX1, InitialX2, InitialX3;
  int BufferIndexH=0, BufferIndexI=0,BufferLength=N_BUFF;
  Xlmk[0]=Xlmk[1]=Xlmk[2]=Xlmk[3]=0.0;
  int nMaxOld = nMax;
  Real epsn = EPS_N;

  // forward loop
  int n=nMax; //start at location of last maximum
  int ContinueSumH=DoSumH, ContinueSumI=DoSumI;
  while(ContinueSumH || ContinueSumI){

    // get new terms
    GetXlmkn(l,m,k,n,ContinueSumH,ContinueSumI);

    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmkn[0]) > TempMaxEIlmkn){
      TempMaxEIlmkn = fabs(Xlmkn[0]);
      nMax = n;
    }
    if(fabs(Xlmkn[0]) > MaxEIlmkn) MaxEIlmkn = fabs(Xlmkn[0]);
    if(fabs(Xlmkn[1]) > MaxEHlmkn) MaxEHlmkn = fabs(Xlmkn[1]);
    if(fabs(Xlmkn[2]) > MaxLIlmkn) MaxLIlmkn = fabs(Xlmkn[2]);
    if(fabs(Xlmkn[3]) > MaxLHlmkn) MaxLHlmkn = fabs(Xlmkn[3]);

    // increment n-index
    n++;

    // check truncation criteria
    ContinueSumH = ( fabs(Xlmkn[1]) > epsn*MaxEHlmkn ||
                     fabs(Xlmkn[3]) > epsn*MaxLHlmkn   );
    ContinueSumI = ( fabs(Xlmkn[0]) > epsn*MaxEIlmkn ||
                     fabs(Xlmkn[2]) > epsn*MaxLIlmkn   );

    // buffer the truncation (H)
    if(DoSumH && BufferIndexH <= BufferLength){
      if(ContinueSumH){
        BufferIndexH=0;
      } else {
        if(BufferIndexH==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1 = fabs(Xlmkn[1]); // EH
          InitialX3 = ( m==0 ? 1.0 : fabs(Xlmkn[3]) ); // LH
        }
        BufferIndexH++;
        ContinueSumH=1;
        if(BufferIndexH > BufferLength){
          ContinueSumH = 0;
          if((InitialX1 < fabs(Xlmkn[1]) ||  // EH modes are increasing
              InitialX3 < fabs(Xlmkn[3]) )   // LH modes are increasing
            && USEnSlope  ) {
            ContinueSumH = 1;
            BufferIndexH = 0;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(DoSumI && BufferIndexI <= BufferLength){
      if(ContinueSumI){
        BufferIndexI=0;
      } else {
        if(BufferIndexI==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0 = fabs(Xlmkn[0]); // EI
          InitialX2 = ( m==0 ? 1.0 : fabs(Xlmkn[2]) ); // LI
        }
        BufferIndexI++;
        ContinueSumI=1;
        if(BufferIndexI > BufferLength){
          ContinueSumI = 0;
          if((InitialX0 < fabs(Xlmkn[0]) ||  // EI modes are increasing
              InitialX2 < fabs(Xlmkn[2]) )  // LI modes are increasing
            && USEnSlope  ) {
            ContinueSumI = 1;
            BufferIndexI = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(n-nMax) > MAXSTEP) Die("MAXSTEP exceeded!");

    // increment output terms
    Xlmk[0] += Xlmkn[0];
    Xlmk[1] += Xlmkn[1];
    Xlmk[2] += Xlmkn[2];
    Xlmk[3] += Xlmkn[3];

  }

  //
  // Because of the symmetry: Z_{l,-m,-k,-n} = (-1)^{l+k} conjugate(Z_{lmkn}),
  // we need only consider positive n.  Since the searches are centered on the 
  // maximum, we usually need to do a backward loop so long as  n stays 
  // non-negative.  
  // 
  n=nMaxOld-1; //start at location of last maximum
  ContinueSumH=DoSumH, ContinueSumI=DoSumI; // make sure the loop starts (for non-negative n that is)
  BufferIndexH = 0, BufferIndexI = 0;
  while((ContinueSumH || ContinueSumI) && n >= 0){

    // get new terms
    GetXlmkn(l,m,k,n,ContinueSumH,ContinueSumI);

    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmkn[0]) > TempMaxEIlmkn){
      TempMaxEIlmkn = fabs(Xlmkn[0]);
      nMax = n;
    }
    if(fabs(Xlmkn[0]) > MaxEIlmkn) MaxEIlmkn = fabs(Xlmkn[0]);
    if(fabs(Xlmkn[1]) > MaxEHlmkn) MaxEHlmkn = fabs(Xlmkn[1]);
    if(fabs(Xlmkn[2]) > MaxLIlmkn) MaxLIlmkn = fabs(Xlmkn[2]);
    if(fabs(Xlmkn[3]) > MaxLHlmkn) MaxLHlmkn = fabs(Xlmkn[3]);

    // increment n-index
    n--;

    // check truncation criteria
    ContinueSumH = ( fabs(Xlmkn[1]) > epsn*MaxEHlmkn ||
                     fabs(Xlmkn[3]) > epsn*MaxLHlmkn   );
    ContinueSumI = ( fabs(Xlmkn[0]) > epsn*MaxEIlmkn ||
                     fabs(Xlmkn[2]) > epsn*MaxLIlmkn  );

    // buffer the truncation (H)
    if(DoSumH && BufferIndexH <= BufferLength){
      if(ContinueSumH){
        BufferIndexH=0;
      } else {
        if(BufferIndexH==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1 = fabs(Xlmkn[1]); // EH
          InitialX3 = ( m==0 ? 1.0 : fabs(Xlmkn[3]) ); // LH
        }
        BufferIndexH++;
        ContinueSumH=1;
        if(BufferIndexH > BufferLength){
          ContinueSumH = 0;
          if((InitialX1 < fabs(Xlmkn[1]) ||  // EH modes are increasing
              InitialX3 < fabs(Xlmkn[3]) )   // LH modes are increasing
            && USEnSlope  ) {
            ContinueSumH = 1;
            BufferIndexH = 0;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(DoSumI && BufferIndexI <= BufferLength){
      if(ContinueSumI){
        BufferIndexI=0;
      } else {
        if(BufferIndexI==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0 = fabs(Xlmkn[0]); // EI
          InitialX2 = ( m==0 ? 1.0 : fabs(Xlmkn[2]) ); // LI
        }
        BufferIndexI++;
        ContinueSumI=1;
        if(BufferIndexI > BufferLength){
          ContinueSumI = 0;
          if((InitialX0 < fabs(Xlmkn[0]) ||  // EI modes are increasing
              InitialX2 < fabs(Xlmkn[2]) )   // LI modes are increasing
            && USEnSlope  ) {
            ContinueSumI = 1;
            BufferIndexI = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(n-nMax) > MAXSTEP) Die("MAXSTEP exceeded!");

    // increment output terms
    Xlmk[0] += Xlmkn[0];
    Xlmk[1] += Xlmkn[1];
    Xlmk[2] += Xlmkn[2];
    Xlmk[3] += Xlmkn[3];

  }

}

void IETD::GetXlmkn(int l, int m, int k, int n, int DoSumH, int DoSumI)
{

  // if this is a circular orbit, only n=0 terms are nonzero.
  if(iekg->e < 1e-2 && n != 0) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  }

  // vanishing frequencies [i.e. (m,k,n) = (0,0,0)] are degenerate cases
  if(fabs(iekg->Omega(m,k,n)) <= 1e-10) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  }

  // TEST ignore the polar voice (enforced elsewhere if theta_ = pi/2) 
  /*
  if(k != 0) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  }
  */

  // TEST ignore the radial voice (enforced elsewhere if e < 0.01)
  /*
  if(n != 0) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  }
  */

  // compute radiation 
  int DoZH = DoSumI;
  int DoZI = DoSumH;
  IEKR iekr(iekg,l,m,k,n,eps,MaxEIlmkn,MaxEHlmkn,DoZH,DoZI);
  Complex ZH = iekr.ZH;
  Complex ZInf = iekr.ZInf;
  Real alpha = iekr.alpha;
  Real omega = iekg->Omega(m,k,n);

  // compute spheroidal harmonic (at a few different angles)
  SWSH swsh(l, -2, m, (iekg->a)*omega);
  Real ct90 = cos(M_PI/2.0 - TEPS);
  Real ct60 = cos(M_PI/3.0);
  Real ct30 = cos(M_PI/6.0);
  Real ct0 = cos(TEPS);
  Real S90 = swsh.spheroid(ct90);
  Real S60 = swsh.spheroid(ct60);
  Real S30 = swsh.spheroid(ct30);
  Real S0 = swsh.spheroid(ct0);

  // compute EInf
  Real EInf = ZH.real()*ZH.real() + ZH.imag()*ZH.imag();
  EInf /= 4.0*M_PI*omega*omega;

  // compute EH
  Real EH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
  EH *= alpha;
  EH /= 4.0*M_PI*omega*omega;

  // compute LInf
  Real LInf = ZH.real()*ZH.real() + ZH.imag()*ZH.imag();
  LInf *= Real(m);
  LInf /= 4.0*M_PI*omega*omega*omega;

  // compute LH
  Real LH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
  LH *= Real(m);
  LH *= alpha;
  LH /= 4.0*M_PI*omega*omega*omega;

  // compute QInf
  Real QInf = (m * iekg->Lz * iekg->avgcot2) + (k * iekg->UpsilonTheta);
  QInf /= -omega;
  QInf += iekg->a * iekg->a * iekg->E * iekg->avgcos2;
  QInf *= -2.0*EInf;

  // compute QH
  Real QH = (m * iekg->Lz * iekg->avgcot2) + (k * iekg->UpsilonTheta);
  QH /= -omega;
  QH += iekg->a * iekg->a * iekg->E * iekg->avgcos2;
  QH *= -2.0*EH;

  // fill output array, remembering that the 
  // n>0 terms should be counted twice.
  Xlmkn[0] = (n > 0 ? 2.0*EInf : EInf);
  Xlmkn[1] = (n > 0 ? 2.0*EH   : EH  );
  Xlmkn[2] = (n > 0 ? 2.0*LInf : LInf);
  Xlmkn[3] = (n > 0 ? 2.0*LH   : LH);

  // record radiation data
  dEH += EH;
  dEI += EInf;
  dLH += LH;
  dLI += LInf;
  //dQ  += (2.0*iekg->Q / iekg->Lz) * (LH + LInf);
  dQ += QInf + QH;
  outfile.precision(16);
  outfile << l << "\t"
       << m << "\t"
       << k << "\t"
       << n << "\t"
       << omega << "\t"
       << ZH.real() << "\t"
       << ZH.imag() << "\t"
       << ZInf.real() <<  "\t"
       << ZInf.imag() << "\t"
       << iekr.alpha << "\t"
       << S90 << "\t"
       << S60 << "\t"
       << S30 << "\t"
       << S0 << "\t"
       << EInf << "\t"
       << EH << "\t"
       << LInf << "\t"
       << LH << "\t"
       //<< QInf << "\t" // TEST
       //<< QH << "\t"  // TEST
       << endl;
  outfile.flush();

  // record the n < 0 entry
  if(n > 0) {
    dEH += EH;
    dEI += EInf;
    dLH += LH;
    dLI += LInf;
    //dQ  += (2.0*iekg->Q / iekg->Lz) * (LH + LInf);
    dQ  += QInf;
    SWSH swsh2(l, -2, -m, -(iekg->a)*omega);
    S90 = swsh2.spheroid(ct90);
    S60 = swsh2.spheroid(ct60);
    S30 = swsh2.spheroid(ct30);
    Real SymFactor = pow(-1.0,Real(l+k));
    outfile << l << "\t"
         << -m << "\t"
         << -k << "\t"
         << -n << "\t"
         << -omega << "\t"
         <<  SymFactor*ZH.real() << "\t"
         << -SymFactor*ZH.imag() << "\t"
         <<  SymFactor*ZInf.real() <<  "\t"
         << -SymFactor*ZInf.imag() << "\t"
         << iekr.alpha << "\t"
         << S90 << "\t"
         << S60 << "\t"
         << S30 << "\t"
         << S0 << "\t"
         << EInf << "\t"
         << EH << "\t"
         << LInf << "\t"
         << LH << "\t"
         // << QInf << "\t" // TEST
         // << QH << "\t"   // TEST
         << endl;
    outfile.flush();
  }

}
