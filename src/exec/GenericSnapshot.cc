// Computes waveforms and fluxes for generic orbits.
// 
// Steve Drasco, August 2007
// 

// standard includes
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <fstream.h>

// Gremlin includes
#include "SNT.h"
#include "SWSH.h"
#include "Globals.h"
#include "IEKG.h"
#include "IEKR.h"

// global varriables.  
//
// read but do not write to the geodesic object
IEKG *iekg; // [a e p theta_] 
// read and write to these during checkpointing
int l_index = 2;
int m_index;
int k_index;
int n_index;
int kMaxOld;
int kMax;
int nMaxOld;
int nMax;
int MAXSTEP = 200;
bool USEnSlope = true;
bool USEkSlope = true;
Real TEPS = 0.001;
Real eps_flux;
Real eps_l;
Real eps_k;
Real eps_n;
Real MaxEIlmkn = 0.0;
Real MaxEIlmk  = 0.0;
Real MaxEIlm   = 0.0;
Real MaxEIl    = 0.0;
Real MaxEHlmkn = 0.0;
Real MaxEHlmk  = 0.0;
Real MaxEHlm   = 0.0;
Real MaxEHl    = 0.0;
Real MaxLIlmkn = 0.0;
Real MaxLIlmk  = 0.0;
Real MaxLIlm   = 0.0;
Real MaxLIl    = 0.0;
Real MaxLHlmkn = 0.0;
Real MaxLHlmk  = 0.0;
Real MaxLHlm   = 0.0;
Real MaxLHl    = 0.0;
Real TempMaxEIlmkn = 0.0;
Real *Xl;    // has elements [0-3]
Real *Xlm;   // has elements [0-3]
Real *Xlmk;  // has elements [0-3]
Real *Xlmkn; // has elements [0-3]
bool ContinueSumH_main;//    = true; 
bool ContinueSumI_main;//    = true;
bool ContinueSumH_GetXl;//   = true;
bool ContinueSumI_GetXl;//   = true;
bool ContinueSumH_GetXlm;//  = true;
bool ContinueSumI_GetXlm;//  = true;
bool ContinueSumH_GetXlmk;// = true;
bool ContinueSumI_GetXlmk;// = true;
int BufferIndexH_GetXlm;
int BufferIndexI_GetXlm;
int BufferLength_GetXlm  = 3;
int BufferIndexH_GetXlmk;
int BufferIndexI_GetXlmk;
int BufferLength_GetXlmk = 5;
Real InitialX0_GetXlm;
Real InitialX1_GetXlm;
Real InitialX2_GetXlm;
Real InitialX3_GetXlm;
Real InitialX0_GetXlmk;
Real InitialX1_GetXlmk;
Real InitialX2_GetXlmk;
Real InitialX3_GetXlmk;
bool RunComplete = false;
// these last few are also global, 
// but are not touched by checkpoint read/write
char fname[500];
fstream outfile;
fstream ckptfile;
bool JUMPSTARTING = false;

// read the checkpoint file
void CheckpointRead(Real a, Real e, Real p, Real MinTheta)
{
  ckptfile.open(fname, ios::in);
  JUMPSTARTING = ckptfile.is_open();
  if(JUMPSTARTING){
    Real theta_ = M_PI * fabs(MinTheta)/180.0;
    Real check_a, check_e, check_p, check_theta_;
    ckptfile >> l_index
             >> m_index
             >> k_index
             >> n_index
             >> kMaxOld
             >> kMax
             >> nMaxOld
             >> nMax
             >> check_a
             >> check_e
             >> check_p
             >> check_theta_
             >> MAXSTEP 
             >> USEnSlope 
             >> USEkSlope
             >> TEPS
             >> eps_flux
             >> eps_l
             >> eps_k
             >> eps_n
             >> MaxEIlmkn
             >> MaxEIlmk
             >> MaxEIlm
             >> MaxEIl 
             >> MaxEHlmkn 
             >> MaxEHlmk
             >> MaxEHlm 
             >> MaxEHl
             >> MaxLIlmkn
             >> MaxLIlmk
             >> MaxLIlm 
             >> MaxLIl 
             >> MaxLHlmkn 
             >> MaxLHlmk 
             >> MaxLHlm 
             >> MaxLHl  
             >> TempMaxEIlmkn
             >> Xl[0]
             >> Xl[1]    
             >> Xl[2]    
             >> Xl[3]    
             >> Xlm[0]    
             >> Xlm[1]
             >> Xlm[2]
             >> Xlm[3]
             >> Xlmk[0]    
             >> Xlmk[1]
             >> Xlmk[2]
             >> Xlmk[3]
             >> Xlmkn[0]
             >> Xlmkn[1]
             >> Xlmkn[2]
             >> Xlmkn[3]
             >> ContinueSumH_main
             >> ContinueSumI_main
             >> ContinueSumH_GetXl
             >> ContinueSumI_GetXl
             >> ContinueSumH_GetXlm
             >> ContinueSumI_GetXlm
             >> ContinueSumH_GetXlmk
             >> ContinueSumI_GetXlmk
             >> BufferIndexH_GetXlm
             >> BufferIndexI_GetXlm
             >> BufferLength_GetXlm
             >> BufferIndexH_GetXlmk
             >> BufferIndexI_GetXlmk
             >> BufferLength_GetXlmk
             >> InitialX0_GetXlm
             >> InitialX1_GetXlm
             >> InitialX2_GetXlm
             >> InitialX3_GetXlm
             >> InitialX0_GetXlmk
             >> InitialX1_GetXlmk
             >> InitialX2_GetXlmk
             >> InitialX3_GetXlmk
             >> RunComplete;
    ckptfile.close();

    // Check the kill flag to see if the run had already finished.
    if(RunComplete){
      outfile.close();
      Die("Checkpoint file indicates that this run is complete. This is a normal exit.");
    }

    bool SameOrbit = ( fabs(a - check_a) / fabs(a)                < 1e-11 ||  
                       fabs(e - check_e) / fabs(e)                < 1e-11 ||  
                       fabs(p - check_p) / fabs(p)                < 1e-11 ||  
                       fabs(theta_ - check_theta_) / fabs(theta_) < 1e-11 );
    if(!SameOrbit){
      Die("The orbital parameters of the input file are not the same as the command line's");
    }
  }
  return;
}

// write to the checkpoint file
void CheckpointWrite()
{
    // open the checkpoint file
    ckptfile.open(fname, ios::out);
    // Quit if we failed to open the file.
    if(!ckptfile.is_open()){
      Die("Failed to open the checkpoint file.");
    }

    // On the JPL supercomputer cosmos, sometimes the error control state
    // of the fstream object ckptfile gets set to "bad", and the program 
    // will not write to the file.  I've found that forcefully resetting the 
    // error control state with clear() causes writing to work as desired.  
    // This isn't a very good solution, but it's the best I have for now.
    ckptfile.clear();

    // write checkpoint data
    ckptfile.precision(16);
    ckptfile       << l_index
           << endl << m_index
           << endl << k_index
           << endl << n_index
           << endl << kMaxOld
           << endl << kMax
           << endl << nMaxOld
           << endl << nMax
           << endl << iekg->a
           << endl << iekg->e
           << endl << iekg->p
           << endl << iekg->theta_
           << endl << MAXSTEP
           << endl << USEnSlope
           << endl << USEkSlope
           << endl << TEPS
           << endl << eps_flux
           << endl << eps_l
           << endl << eps_k
           << endl << eps_n
           << endl << MaxEIlmkn
           << endl << MaxEIlmk
           << endl << MaxEIlm
           << endl << MaxEIl
           << endl << MaxEHlmkn
           << endl << MaxEHlmk
           << endl << MaxEHlm
           << endl << MaxEHl
           << endl << MaxLIlmkn
           << endl << MaxLIlmk
           << endl << MaxLIlm
           << endl << MaxLIl
           << endl << MaxLHlmkn
           << endl << MaxLHlmk
           << endl << MaxLHlm
           << endl << MaxLHl
           << endl << TempMaxEIlmkn
           << endl << Xl[0]    
           << endl << Xl[1]
           << endl << Xl[2]
           << endl << Xl[3]
           << endl << Xlm[0]
           << endl << Xlm[1]
           << endl << Xlm[2]
           << endl << Xlm[3]
           << endl << Xlmk[0]
           << endl << Xlmk[1]
           << endl << Xlmk[2]
           << endl << Xlmk[3]
           << endl << Xlmkn[0]
           << endl << Xlmkn[1]
           << endl << Xlmkn[2]
           << endl << Xlmkn[3]
           << endl << ContinueSumH_main
           << endl << ContinueSumI_main
           << endl << ContinueSumH_GetXl
           << endl << ContinueSumI_GetXl
           << endl << ContinueSumH_GetXlm
           << endl << ContinueSumI_GetXlm
           << endl << ContinueSumH_GetXlmk
           << endl << ContinueSumI_GetXlmk
           << endl << BufferIndexH_GetXlm
           << endl << BufferIndexI_GetXlm
           << endl << BufferLength_GetXlm
           << endl << BufferIndexH_GetXlmk
           << endl << BufferIndexI_GetXlmk
           << endl << BufferLength_GetXlmk
           << endl << InitialX0_GetXlm
           << endl << InitialX1_GetXlm
           << endl << InitialX2_GetXlm
           << endl << InitialX3_GetXlm
           << endl << InitialX0_GetXlmk
           << endl << InitialX1_GetXlmk
           << endl << InitialX2_GetXlmk
           << endl << InitialX3_GetXlmk
           << endl << RunComplete
           << endl;
  ckptfile.flush();
  ckptfile.close();
  return;
}

// subroutines for main loop
void GetXlmkn()
{

  // if this is a circular orbit, only n=0 terms are nonzero.
  if(iekg->e < 1e-2 && n_index != 0) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  }

  // vanishing frequencies [i.e. (m,k,n) = (0,0,0)] are degenerate cases
  if(fabs(iekg->Omega(m_index,k_index,n_index)) <= 1e-10) {
    Xlmkn[0]=Xlmkn[1]=Xlmkn[2]=Xlmkn[3]=0.0;
    return;
  } 

  // If we were jumpstarting a previous calculation, then the jumpstarting is now done.
  if(JUMPSTARTING) {
    JUMPSTARTING = false;
    // We return with the Xlmkn that were computed before, and so that the 
    // preceeding routines can decide how to increment the indices.
    return;
  }

  //compute radiation 
  int DoZH = ContinueSumI_GetXlmk;
  int DoZI = ContinueSumH_GetXlmk;
  IEKR iekr(iekg,l_index,m_index,k_index,n_index,eps_flux,MaxEIlmkn,MaxEHlmkn,DoZH,DoZI);
  Complex ZH = iekr.ZH;
  Complex ZInf = iekr.ZInf;
  Complex Ain = iekr.Ain;
  Real alpha = iekr.alpha;
  Real omega = iekg->Omega(m_index,k_index,n_index);

  // compute spheroidal harmonic (at a few different angles)
  SWSH swsh(l_index, -2, m_index, (iekg->a)*omega);
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
  LInf *= Real(m_index);
  LInf /= 4.0*M_PI*omega*omega*omega;
                                                                                                   
  // compute LH
  Real LH = ZInf.real()*ZInf.real() + ZInf.imag()*ZInf.imag();
  LH *= Real(m_index);
  LH *= alpha;
  LH /= 4.0*M_PI*omega*omega*omega;

  // compute QInf
  Real QInf = (m_index * iekg->Lz * iekg->avgcot2) + (k_index * iekg->UpsilonTheta);
  QInf /= omega;
  QInf -= iekg->a * iekg->a * iekg->E * iekg->avgcos2;
  QInf *= 2.0*EInf;

  // compute QH
  Real QH = (m_index * iekg->Lz * iekg->avgcot2) + (k_index * iekg->UpsilonTheta);
  QH /= omega;
  QH -= iekg->a * iekg->a * iekg->E * iekg->avgcos2;
  QH *= 2.0*EH;

  // fill output array, remembering that the 
  // n>0 terms should be counted twice.
  Xlmkn[0] = (n_index > 0 ? 2.0*EInf : EInf);
  Xlmkn[1] = (n_index > 0 ? 2.0*EH   : EH  );
  Xlmkn[2] = (n_index > 0 ? 2.0*LInf : LInf);
  Xlmkn[3] = (n_index > 0 ? 2.0*LH   : LH);
  
  // record radiation data
  outfile.precision(16); 
  outfile << l_index << "\t"
       << m_index << "\t"
       << k_index << "\t"
       << n_index << "\t"
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
       //<< QH << "\t"   // TEST
       << Ain.real() << "\t"
       << Ain.imag() << "\t"
       << endl;
  outfile.flush();

  // record the n < 0 entry
  if(n_index > 0) {
    SWSH swsh2(l_index, -2, -m_index, -(iekg->a)*omega);
    S90 = swsh2.spheroid(ct90);
    S60 = swsh2.spheroid(ct60);
    S30 = swsh2.spheroid(ct30);
    Real SymFactor = pow(-1.0,Real(l_index+k_index));
    outfile << l_index << "\t"
         << -m_index << "\t"
         << -k_index << "\t"
         << -n_index << "\t"
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
         //<< QInf << "\t" // TEST
         //<< QH << "\t"   // TEST
         << Ain.real() << "\t"
         << Ain.imag() << "\t"
         << endl;
    outfile.flush();
  } 

  // TEST
  /*
  //outfile << "this is output sent from GenericSnapshot to the output file.\n";
  //cout << "this is output sent from GenericSnapshot to cout.\n";
  ckptfile.open(fname, ios::out);
  ckptfile << "this is output sent from GenericSnapshot to the checkpoint file.\n";
  cout << "ckptfile.rdstate() = " << ckptfile.rdstate() << endl;
  cout << "ckptfile.good() = " << ckptfile.good() << endl;
  ckptfile.clear();
  ckptfile << "this is output sent from GenericSnapshot to the checkpoint file.\n";
  cout << "ckptfile.rdstate() = " << ckptfile.rdstate() << endl;
  cout << "ckptfile.good() = " << ckptfile.good() << endl;
  ckptfile.close();
  Die("Test Die");
  */
  // END TEST

  // Update the checkpoint file 
  CheckpointWrite(); 

}
void GetXlmk(){ 

  // if this is an equatorial orbit, only k=0 terms are nonzero. 
  if(iekg->theta_ == 0.5*M_PI && k_index != 0) {
    Xlmk[0]=Xlmk[1]=Xlmk[2]=Xlmk[3]=0.0;
    return;
  }

  // Begin the forward loop.  If jumpstarting, only do this loop if the 
  // calculation we are jumpstarting was actually in this loop when it was 
  // interupted. 
  if( !JUMPSTARTING || (JUMPSTARTING && n_index >= nMaxOld) ) {

    // initialize the sum over k, unless we are jumpstarting
    if(!JUMPSTARTING){
      Xlmk[0]=Xlmk[1]=Xlmk[2]=Xlmk[3]=0.0;
      n_index = nMax;
      nMaxOld = nMax;
      ContinueSumH_GetXlmk = ContinueSumH_GetXlm;
      ContinueSumI_GetXlmk = ContinueSumI_GetXlm;
      BufferIndexH_GetXlmk = 0;
      BufferIndexI_GetXlmk = 0;
    }
    while(ContinueSumH_GetXlmk || ContinueSumI_GetXlmk){
  
      // get new terms
      GetXlmkn();
  
      // check for any new global maxima
      // we track the peak for EI
      if(fabs(Xlmkn[0]) > TempMaxEIlmkn){
        TempMaxEIlmkn = fabs(Xlmkn[0]);
        nMax = n_index;
      }
      if(fabs(Xlmkn[0]) > MaxEIlmkn) MaxEIlmkn = fabs(Xlmkn[0]);
      if(fabs(Xlmkn[1]) > MaxEHlmkn) MaxEHlmkn = fabs(Xlmkn[1]);
      if(fabs(Xlmkn[2]) > MaxLIlmkn) MaxLIlmkn = fabs(Xlmkn[2]);
      if(fabs(Xlmkn[3]) > MaxLHlmkn) MaxLHlmkn = fabs(Xlmkn[3]);

      // increment n-index
      n_index++;
                                                                                                                             
      // check truncation criteria
      ContinueSumH_GetXlmk = ( fabs(Xlmkn[1]) > eps_n*MaxEHlmkn ||
                               fabs(Xlmkn[3]) > eps_n*MaxLHlmkn   );
      ContinueSumI_GetXlmk = ( fabs(Xlmkn[0]) > eps_n*MaxEIlmkn ||
                               fabs(Xlmkn[2]) > eps_n*MaxLIlmkn   );

      // buffer the truncation (H)
      if(ContinueSumH_GetXlm && BufferIndexH_GetXlmk <= BufferLength_GetXlmk){
        if(ContinueSumH_GetXlmk){
          BufferIndexH_GetXlmk=0;
        } else {
          if(BufferIndexH_GetXlmk==0){
            // This is the start of the buffer.  We set initial values, that will
            // be used to make sure modes are decreasing when we truncate.
            InitialX1_GetXlmk = fabs(Xlmkn[1]); // EH
            InitialX3_GetXlmk = ( m_index==0 ? 1.0 : fabs(Xlmkn[3]) ); // LH
          }
          BufferIndexH_GetXlmk++;
          ContinueSumH_GetXlmk=true;
          if(BufferIndexH_GetXlmk > BufferLength_GetXlmk){
            ContinueSumH_GetXlmk = false;
            if((InitialX1_GetXlmk < fabs(Xlmkn[1]) ||  // EH modes are increasing
                InitialX3_GetXlmk < fabs(Xlmkn[3]) )   // LH modes are increasing
              && USEnSlope  ) {
              ContinueSumH_GetXlmk = true;
              BufferIndexH_GetXlmk = 0;
            }
          }
        }
      }

      // buffer the truncation (I)
      if(ContinueSumI_GetXlm && BufferIndexI_GetXlmk <= BufferLength_GetXlmk){
        if(ContinueSumI_GetXlmk){
          BufferIndexI_GetXlmk=0;
        } else {
          if(BufferIndexI_GetXlmk==0){
            // This is the start of the buffer.  We set initial values, that will
            // be used to make sure modes are decreasing when we truncate.
            InitialX0_GetXlmk = fabs(Xlmkn[0]); // EI
            InitialX2_GetXlmk = ( m_index==0 ? 1.0 : fabs(Xlmkn[2]) ); // LI
          }
          BufferIndexI_GetXlmk++;
          ContinueSumI_GetXlmk=true;
          if(BufferIndexI_GetXlmk > BufferLength_GetXlmk){
            ContinueSumI_GetXlmk = false;
            if((InitialX0_GetXlmk < fabs(Xlmkn[0]) ||  // EI modes are increasing
                InitialX2_GetXlmk < fabs(Xlmkn[2]) )  // LI modes are increasing
              && USEnSlope  ) {
              ContinueSumI_GetXlmk = true;
              BufferIndexI_GetXlmk = 0;
            }
          }
        }
      }
 
      // we only take so many steps
      if(abs(n_index-nMax) > MAXSTEP) Die("MAXSTEP exceeded!");

      // increment output terms
      Xlmk[0] += Xlmkn[0];
      Xlmk[1] += Xlmkn[1];
      Xlmk[2] += Xlmkn[2];
      Xlmk[3] += Xlmkn[3];

    }

  } // end positive loop 

  //
  // Because of the symmetry: Z_{l,-m,-k,-n} = (-1)^{l+k} conjugate(Z_{lmkn}),
  // we need only consider positive n.  Since the searches are centered on the 
  // maximum, we usually need to do a backward loop so long as  n stays 
  // non-negative.  
  // 
  // Begin or resume the backward loop.  If we are not jumpstarting, 
  // some things need to be initialized.
  if(!JUMPSTARTING){
    n_index=nMaxOld-1; //start at location of last maximum
    ContinueSumH_GetXlmk = ContinueSumH_GetXlm;
    ContinueSumI_GetXlmk = ContinueSumI_GetXlm; // make sure the loop starts (for non-negative n that is)
    BufferIndexH_GetXlmk = 0;
    BufferIndexI_GetXlmk = 0;
  }
  while((ContinueSumH_GetXlmk || ContinueSumI_GetXlmk) && n_index >= 0){

    // get new terms
    GetXlmkn();
    
    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmkn[0]) > TempMaxEIlmkn){
      TempMaxEIlmkn = fabs(Xlmkn[0]);
      nMax = n_index;
    }
    if(fabs(Xlmkn[0]) > MaxEIlmkn) MaxEIlmkn = fabs(Xlmkn[0]);
    if(fabs(Xlmkn[1]) > MaxEHlmkn) MaxEHlmkn = fabs(Xlmkn[1]);
    if(fabs(Xlmkn[2]) > MaxLIlmkn) MaxLIlmkn = fabs(Xlmkn[2]);
    if(fabs(Xlmkn[3]) > MaxLHlmkn) MaxLHlmkn = fabs(Xlmkn[3]);
                                                                                                                                                             
    // increment n-index
    n_index--;

    // check truncation criteria
    ContinueSumH_GetXlmk = ( fabs(Xlmkn[1]) > eps_n*MaxEHlmkn ||
                             fabs(Xlmkn[3]) > eps_n*MaxLHlmkn   );
    ContinueSumI_GetXlmk = ( fabs(Xlmkn[0]) > eps_n*MaxEIlmkn ||
                             fabs(Xlmkn[2]) > eps_n*MaxLIlmkn  );

    // buffer the truncation (H)
    if(ContinueSumH_GetXlm && BufferIndexH_GetXlmk <= BufferLength_GetXlmk){
      if(ContinueSumH_GetXlmk){
        BufferIndexH_GetXlmk=0;
      } else {
        if(BufferIndexH_GetXlmk==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1_GetXlmk = fabs(Xlmkn[1]); // EH
          InitialX3_GetXlmk = ( m_index==0 ? 1.0 : fabs(Xlmkn[3]) ); // LH
        }
        BufferIndexH_GetXlmk++;
        ContinueSumH_GetXlmk=true;
        if(BufferIndexH_GetXlmk > BufferLength_GetXlmk){
          ContinueSumH_GetXlmk = false;
          if((InitialX1_GetXlmk < fabs(Xlmkn[1]) ||  // EH modes are increasing
              InitialX3_GetXlmk < fabs(Xlmkn[3]) )   // LH modes are increasing
            && USEnSlope  ) {
            ContinueSumH_GetXlmk = true;
            BufferIndexH_GetXlmk = false;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(ContinueSumI_GetXlm && BufferIndexI_GetXlmk <= BufferLength_GetXlmk){
      if(ContinueSumI_GetXlmk){
        BufferIndexI_GetXlmk=0;
      } else {
        if(BufferIndexI_GetXlmk==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0_GetXlmk = fabs(Xlmkn[0]); // EI
          InitialX2_GetXlmk = ( m_index==0 ? 1.0 : fabs(Xlmkn[2]) ); // LI
        }
        BufferIndexI_GetXlmk++;
        ContinueSumI_GetXlmk=true;
        if(BufferIndexI_GetXlmk > BufferLength_GetXlmk){
          ContinueSumI_GetXlmk = false;
          if((InitialX0_GetXlmk < fabs(Xlmkn[0]) ||  // EI modes are increasing
              InitialX2_GetXlmk < fabs(Xlmkn[2]) )   // LI modes are increasing
            && USEnSlope  ) {
            ContinueSumI_GetXlmk = true;
            BufferIndexI_GetXlmk = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(n_index-nMax) > MAXSTEP) Die("MAXSTEP exceeded!");

    // increment output terms
    Xlmk[0] += Xlmkn[0];
    Xlmk[1] += Xlmkn[1];
    Xlmk[2] += Xlmkn[2];
    Xlmk[3] += Xlmkn[3];

  }

}
void GetXlm(){ 

  // Begin the forward loop.  If jumpstarting, only do this loop if the 
  // calculation we are jumpstarting was actually in this loop when it was 
  // interupted. 
  if( !JUMPSTARTING || (JUMPSTARTING && k_index >= kMaxOld) ) {

    // initialize the sum over k, unless we are jumpstarting
    if(!JUMPSTARTING){
      Xlm[0]=Xlm[1]=Xlm[2]=Xlm[3]=0.0;
      k_index = kMax;
      kMaxOld = kMax;
      ContinueSumH_GetXlm = ContinueSumH_GetXl;
      ContinueSumI_GetXlm = ContinueSumI_GetXl;
      BufferIndexH_GetXlm = 0;
      BufferIndexI_GetXlm = 0;
    }

    while(ContinueSumH_GetXlm || ContinueSumI_GetXlm){

      // get new terms
      GetXlmk();

      // check for any new global maxima
      // we track the peak for EI
      if(fabs(Xlmk[0]) > MaxEIlmk){
        MaxEIlmk = fabs(Xlmk[0]);
        kMax = k_index;
      }
      if(fabs(Xlmk[1]) > MaxEHlmk) MaxEHlmk = fabs(Xlmk[1]);
      if(fabs(Xlmk[2]) > MaxLIlmk) MaxLIlmk = fabs(Xlmk[2]);
      if(fabs(Xlmk[3]) > MaxLHlmk) MaxLHlmk = fabs(Xlmk[3]);

      // increment k-index
      k_index++;

      // check truncation criteria
      ContinueSumH_GetXlm = ( fabs(Xlmk[1]) > eps_k*MaxEHlmk ||
                       fabs(Xlmk[3]) > eps_k*MaxLHlmk   );
      ContinueSumI_GetXlm = ( fabs(Xlmk[0]) > eps_k*MaxEIlmk ||
                       fabs(Xlmk[2]) > eps_k*MaxLIlmk   );

      // buffer the truncation (H)
      if(ContinueSumH_GetXl && BufferIndexH_GetXlm <= BufferLength_GetXlm){
        if(ContinueSumH_GetXlm){
          BufferIndexH_GetXlm=0;
        } else {
          if(BufferIndexH_GetXlm==0){
            // This is the start of the buffer.  We set initial values, that will
            // be used to make sure modes are decreasing when we truncate.
            InitialX1_GetXlm = fabs(Xlmk[1]); // EH
            InitialX3_GetXlm = ( m_index==0 ? 1.0 : fabs(Xlmk[3]) ); // LH
          }
          BufferIndexH_GetXlm++;
          ContinueSumH_GetXlm=ContinueSumH_GetXl; 
          if(BufferIndexH_GetXlm > BufferLength_GetXlm){
            ContinueSumH_GetXlm = false;
            if((InitialX1_GetXlm < fabs(Xlmk[1]) ||  // EH modes are increasing
                InitialX3_GetXlm < fabs(Xlmk[3]) )   // LH modes are increasing
              && USEkSlope  ) {
              ContinueSumH_GetXlm = true;
              BufferIndexH_GetXlm = 0;
            }
          } 
        }
      }

      // buffer the truncation (I)
      if(ContinueSumI_GetXl && BufferIndexI_GetXlm <= BufferLength_GetXlm){
        if(ContinueSumI_GetXlm){
          BufferIndexI_GetXlm=0;
        } else {
          if(BufferIndexI_GetXlm==0){
            // This is the start of the buffer.  We set initial values, that will
            // be used to make sure modes are decreasing when we truncate.
            InitialX0_GetXlm = fabs(Xlmk[0]); // EI
            InitialX2_GetXlm = ( m_index==0 ? 1.0 : fabs(Xlmk[2]) ); // LI
          }
          BufferIndexI_GetXlm++;
          ContinueSumI_GetXlm=ContinueSumI_GetXl;
          if(BufferIndexI_GetXlm > BufferLength_GetXlm){
            ContinueSumI_GetXlm = false;
            if((InitialX0_GetXlm < fabs(Xlmk[0]) ||  // EI modes are increasing
                InitialX2_GetXlm < fabs(Xlmk[2]) )   // LI modes are increasing
              && USEkSlope  ) {
              ContinueSumI_GetXlm = true;
              BufferIndexI_GetXlm = 0;
            }
          }
        }
      }

      // we only take so many steps
      if(abs(k_index-kMax) > MAXSTEP) Die("MAXSTEP exceeded!");
                                                                                                                             
      // increment output terms
      Xlm[0] += Xlmk[0];
      Xlm[1] += Xlmk[1];
      Xlm[2] += Xlmk[2];
      Xlm[3] += Xlmk[3];

    }

  }  // end of forward loop 

  // Begin or resume the backward loop.  If we are not jumpstarting, 
  // some things need to be initialized.
  if(!JUMPSTARTING){
    k_index = kMaxOld-1;  // start at location of last maximum
    ContinueSumH_GetXlm = ContinueSumH_GetXl;
    ContinueSumI_GetXlm = ContinueSumI_GetXl; // make sure the loop starts
    BufferIndexH_GetXlm = 0;
    BufferIndexI_GetXlm = 0;
  }
  while(ContinueSumH_GetXlm || ContinueSumI_GetXlm){

    // get new terms
    GetXlmk();

    // check for any new global maxima
    // we track the peak for EI
    if(fabs(Xlmk[0]) > MaxEIlmk){
      MaxEIlmk = fabs(Xlmk[0]);
      kMax = k_index;
    }
    if(fabs(Xlmk[1]) > MaxEHlmk) MaxEHlmk = fabs(Xlmk[1]);
    if(fabs(Xlmk[2]) > MaxLIlmk) MaxLIlmk = fabs(Xlmk[2]);
    if(fabs(Xlmk[3]) > MaxLHlmk) MaxLHlmk = fabs(Xlmk[3]);
                                                                                                                                                             
    // increment k-index
    k_index--;

    // check truncation criteria
    ContinueSumH_GetXlm = ( fabs(Xlmk[1]) > eps_k*MaxEHlmk ||
                            fabs(Xlmk[3]) > eps_k*MaxLHlmk   );
    ContinueSumI_GetXlm = ( fabs(Xlmk[0]) > eps_k*MaxEIlmk ||
                            fabs(Xlmk[2]) > eps_k*MaxLIlmk   );

    // buffer the truncation (H)
    if(ContinueSumH_GetXl && BufferIndexH_GetXlm <= BufferLength_GetXlm){
      if(ContinueSumH_GetXlm){
        BufferIndexH_GetXlm=0;
      } else {
        if(BufferIndexH_GetXlm==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX1_GetXlm = fabs(Xlmk[1]); // EH
          InitialX3_GetXlm = ( m_index==0 ? 1.0 : fabs(Xlmk[3]) ); // LH
        }
        BufferIndexH_GetXlm++;
        ContinueSumH_GetXlm=true;
        if(BufferIndexH_GetXlm > BufferLength_GetXlm){
          ContinueSumH_GetXlm = false;
          if((InitialX1_GetXlm < fabs(Xlmk[1]) ||  // EH modes are increasing
              InitialX3_GetXlm < fabs(Xlmk[3]) )   // LH modes are increasing
            && USEkSlope  ) {
            ContinueSumH_GetXlm = true;
            BufferIndexH_GetXlm = 0;
          }
        }
      }
    }

    // buffer the truncation (I)
    if(ContinueSumI_GetXl && BufferIndexI_GetXlm <= BufferLength_GetXlm){
      if(ContinueSumI_GetXlm){
        BufferIndexI_GetXlm=0;
      } else {
        if(BufferIndexI_GetXlm==0){
          // This is the start of the buffer.  We set initial values, that will
          // be used to make sure modes are decreasing when we truncate.
          InitialX0_GetXlm = fabs(Xlmk[0]); // EI
          InitialX2_GetXlm = ( m_index==0 ? 1.0 : fabs(Xlmk[2]) ); // LI
        }
        BufferIndexI_GetXlm++;
        ContinueSumI_GetXlm=true;
        if(BufferIndexI_GetXlm > BufferLength_GetXlm){
          ContinueSumI_GetXlm = false;
          if((InitialX0_GetXlm < fabs(Xlmk[0]) ||  // EI modes are increasing
              InitialX2_GetXlm < fabs(Xlmk[2]) )   // LI modes are increasing
            && USEkSlope  ) {
            ContinueSumI_GetXlm = true;
            BufferIndexI_GetXlm = 0;
          }
        }
      }
    }

    // we only take so many steps
    if(abs(k_index-kMax) > MAXSTEP) Die("MAXSTEP exceeded!");
                                                                                                                             
    // increment output terms
    Xlm[0] += Xlmk[0];
    Xlm[1] += Xlmk[1];
    Xlm[2] += Xlmk[2];
    Xlm[3] += Xlmk[3];

  }

}
void GetXl(){

  // initialize the sum over m, unless we are jumpstarting
  if(!JUMPSTARTING){
    Xl[0]=Xl[1]=Xl[2]=Xl[3]=0.0;
    m_index = -l_index;
    // We do not truncate the m-sum.  However, defining
    // "ContinueSumH_GetXl" and "ContinueSumI_GetXl" will clairfy 
    // the subroutines below if truncation of the m-sum is ever added.
    ContinueSumH_GetXl = ContinueSumH_main;
    ContinueSumI_GetXl = ContinueSumI_main;
  }

  //
  // go through all possible values for m
  //
  for(m_index; m_index<=l_index; m_index++){
    
    // get new terms
    GetXlm();

    // increment output terms
    Xl[0] += Xlm[0];
    Xl[1] += Xlm[1];
    Xl[2] += Xlm[2];
    Xl[3] += Xlm[3];
  }

}


int main(int argc, char **argv)
{

  // read input
  ios::sync_with_stdio(); // call this function so iostraem works with stdio
  if(argc != 7) {
    cerr << "Arguments: 1. a  2. e  3. p  4. MinTheta (degrees) 5. eps 6. OutputFilename" << endl;
    exit(1);
  }
  Real a = atof(argv[1]);
  Real e = atof(argv[2]);
  Real p = atof(argv[3]);
  Real MinTheta = atof(argv[4]);
  eps_flux = atof(argv[5]);
  sprintf(fname, "%s", argv[6]);

  // open an output file
  //
  //  This is the main data file where we dump all the information 
  //  that describes the snapshot we are working on.
  //
  //  Use ios::app so that data is only appended to the file, since 
  //  we might only be suplimenting an existing data file from a 
  //  previous calculation that was interupted.
  //
  outfile.open(fname, ios::out | ios::app);
  // Quit if we failed to open the file.
  if(!outfile.is_open()){
    Die("Failed to open the output file.");
  }

  // convert fname from the name of the main output file to the name of the 
  // checkpoint file by adding extension .ckpt
  strcat(fname, ".ckpt");

  // initialize global vectors
  Xl    = Realvector(0,3);
  Xlm   = Realvector(0,3);
  Xlmk  = Realvector(0,3);
  Xlmkn = Realvector(0,3);

  // open a ckeckpoint file (same file name as the main output file name, 
  // but with extension .ckpt)
  //
  //  This is a suplimental data file where we store information 
  //  describing the current status of the calculation.  
  //
  //  If the checkpoint file already exists, it will be read and a calculation 
  //  will be "jumpstarted" !so that it picks up where the last one ended.  
  //  If the file doesn't exist, it will be created so that future runs can 
  //  "jumpstart" the calculation we are about to begin, should it be interupted.
  //
  CheckpointRead(a, e, p, MinTheta);
  

  // Set initial nMax using Peters-Mathews formula.
  // If we are jumpstarting a previous calculation, 
  // don't overwtire the existing nMax.
  if(e > 0.1 & !JUMPSTARTING) {
    nMax = int( exp( 0.5 - 1.5*log(1 - e)) );
  }

  // define the accuracies
  eps_l = eps_flux;
  eps_k = eps_flux;
  eps_n = eps_flux;

  // translate MinTheta to [theta_ (radians) & prograde (integer)]
  Real theta_ = M_PI * fabs(MinTheta)/180.0;
  int prograde = (MinTheta > 0 ? 1 : 0);

  // compute the geodesic
  IEKG iekg_tmp(prograde,a,e,p,theta_,1e-12);
  iekg = &iekg_tmp;
  
  // If this is a new calculation, record geodesic information 
  // as matlab-commented header in the main output file.
  if(!JUMPSTARTING){
    outfile.precision(16);
    outfile << "%\n"
       << "%r_isco/M = " << iekg->r_isco << "\n"
       << "%a = " << a << "\n"
       << "%e = " << e << "\n"
       << "%p = " << p << "\n"
       << "%MinTheta = " << MinTheta << " (degrees)\n"
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
       << "%eps = " << eps_flux << "\n"
       << "%eps_l = " << eps_l << "\n"
       << "%eps_k = " << eps_k << "\n"
       << "%eps_n = " << eps_n << "\n"
       << "%Teps = " << TEPS << "\n"
       << "%Delta_t_kmax = " << iekg->Delta_t_kmax << "\n"
       << "%Delta_t_nmax = " << iekg->Delta_t_nmax << "\n"
       << "%Delta_phi_kmax = " << iekg->Delta_phi_kmax << "\n"
       << "%Delta_phi_nmax = " << iekg->Delta_phi_nmax << "\n"
       << "%MAXSTEP = " << MAXSTEP << "\n"
       << "%K_BUFF = " << BufferLength_GetXlm << "\n"
       << "%N_BUFF = " << BufferLength_GetXlmk << "\n"
       << "%USEnSlope = " << USEnSlope << "\n"
       << "%USEkSlope = " << USEkSlope << "\n"
       << "%\n"
       << "%below: [l m k n omega Z^H.real() ZH.imag() Z^Inf.real() ZInf.imag() ... \n"
       // << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH QInf QH Ain.real() Ain.imag()]" << "\n" // TEST
       << "%    ... alpha S(pi/2-Teps) S(pi/3) S(pi/6) S(Teps) EInf EH LInf LH Ain.real() Ain.imag()]" << "\n"
       << "%\n";
    outfile.flush();
  }

  // if note jumpstarting, initialize the main sum over l
  if(!JUMPSTARTING){
    l_index = 2;
    ContinueSumH_main = false; // false = don't compute horizon fluxes.  true = compute horizon fluxes.
    ContinueSumI_main = true;
  }

  // outermost loop over l-index is done here
  // data is recorded internaly
  while(ContinueSumH_main || ContinueSumI_main){

    // reset temporary peak value of (dE/dt)^{\infty}
    // for tracking peaks in n-index, since the peak's location 
    // looks linear in l-index.  Skip this step if jumpstarting.
    if(!JUMPSTARTING){
      TempMaxEIlmkn = 0.0;
    }

    // get next set of terms in sum
    GetXl();

    // check if we have any new maxima
    if(fabs(Xl[0]) > MaxEIl) MaxEIl = fabs(Xl[0]);
    if(fabs(Xl[1]) > MaxEHl) MaxEHl = fabs(Xl[1]);
    if(fabs(Xl[2]) > MaxLIl) MaxLIl = fabs(Xl[2]);
    if(fabs(Xl[3]) > MaxLHl) MaxLHl = fabs(Xl[3]);

    // increment l-index
    l_index++;

    // check truncation criteria
    ContinueSumH_main = (fabs(Xl[1]) > eps_l*MaxEHl ||
                         fabs(Xl[3]) > eps_l*MaxLHl   );
    ContinueSumI_main = (fabs(Xl[0]) > eps_l*MaxEIl ||
                         fabs(Xl[2]) > eps_l*MaxLIl   );

    // we only take so many steps
    if(l_index-2 > MAXSTEP) Die("MAXSTEP exceeded!");

  }

  // write a closing line to the data file
  outfile << "%\n"  
          << "% normal exit - calculation complete \n";
  outfile.flush();

  // close the main output file
  outfile.close();

  // Tell the checkpoint file that the run is finished so that future 
  // attempts to jumpstart the run will know to quit.
  RunComplete = true;
  CheckpointWrite();

}
