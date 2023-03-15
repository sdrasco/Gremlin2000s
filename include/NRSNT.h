//---------------------------------------------------------------------------
//
// $Id: NRSNT.h,v 1.1.1.1 2003/08/09 19:19:39 sdrasco Exp $
//
//---------------------------------------------------------------------------
//
// Prototypes for NR routines needed for class SNT.
//
#ifndef _NRSNT_H
#define _NRSNT_H

#include "Globals.h"

void pzextr(int iest, Real xest, Complex yest[], Complex yz[],
            Complex dy[], Complex **d, Real *x, int nv);
void rzextr(int iest, Real xest, Complex yest[], Complex yz[],
            Complex dy[], Complex **d, Real *x, int nv);
void rzextr(int iest, Real xest, Real yest[], Real yz[],
            Real dy[], Real **d, Real *x, int nv);

#endif
