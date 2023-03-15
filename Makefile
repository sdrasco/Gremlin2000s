#############################################################################
##
## $Id: Makefile,v 1.47 2008/12/04 00:40:40 sdrasco Exp $
##
#############################################################################
#
# Makefile for Pebble.
#
#############################################################################

.PHONY : dummy

#############################################################################
#
# All the paths we need to know.
#
#############################################################################

TOP = $(HOME)/Proj/Gremlin
INCL = $(TOP)/include
BIN = $(TOP)/bin
LIB = $(TOP)/lib

SRC = $(TOP)/src
CIRCSRC = $(SRC)/circ
CIRCEQSRC = $(SRC)/circeq
EXECSRC = $(SRC)/exec
IESRC = $(SRC)/inclecc
NRSRC = $(SRC)/numrec
SNTSRC = $(SRC)/sasnakteuk
SWSHSRC = $(SRC)/swsh
UTILSRC = $(SRC)/utility
ALLSRCS = $(CIRCSRC):$(CIRCEQSRC):$(EXECSRC):$(IESRC):$(NRSRC):$(SNTSRC):$(SWSHSRC):$(UTILSRC)

VPATH = $(BIN):$(INCL):$(LIB):$(ALLSRCS)

#############################################################################
#
# Useful variables.
#
#############################################################################

CC = g++

AR = ar rv

SYSLIBS = -lm 

# generic
#CFLAGS = -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated
#CFLAGS = -g -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated
#CFLAGS = -g -pg -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated
#CFLAGS = -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated

# intel core duo
#CFLAGS = -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated -march=prescott

# intel core 2 duo
CFLAGS = -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated -march=nocona

# Xeon
#CFLAGS =  -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated -fomit-frame-pointer -malign-double -march=pentium4

# AMD
#CFLAGS =  -O3 -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated -fomit-frame-pointer -malign-double -march=athlon-xp

# for G5
#CFLAGS =  -fast -Wall -Wno-unused -Wno-uninitialized -Wno-deprecated

#############################################################################
#
# -Wall to catch as many warnings as possible, except:
#
# -Wno-unused because some Numerical Recipes functions don't use all of
#             their parameters, and
#
# -Wno-uninitialized because a few Numerical Recipes functions have
#                    variables that look like they'll be used uninitialized.
#                    (They won't be, don't worry about it.)
#
# -Wno-deprecated because on newer compilers, my C++ generates complaints
#                 about the header format.
#
#############################################################################

#############################################################################
#
# All .o files.
#
#############################################################################

CIRCOBJS = CID.o CKG.o CKR.o CTD.o

CIRCEQOBJS = CEID.o CEKG.o CEKR.o CETD.o

IEOBJS = IEKG.o IEKR.o IETD.o

GKGOBJS = GKG.o

INTOBJS = Integral.o

LBOBJS = LB.o

NROBJS =NRElle.o NREllf.o NREllpi.o NRFactrl.o NRGammln.o NRIndexx.o \
	NRLUbksb.o NRLUdcmp.o NRPythag.o NRUtil.o \
	NRpzextr.o NRrc.o NRrd.o NRrf.o NRrj.o \
	NRrzextr.o NRspline.o NRsplint.o NRtqli.o NRtred2.o 

RRGWOBJS = RRGW.o

SNTOBJS = SNTBSstep.o SNTMmid.o SNTOdeint.o SNT.o SNTSasNak.o SNTTeuk.o SNTInterp.o

SWSHOBJS = SWSHCGUtil.o SWSHSpherical.o SWSHSpheroid.o

.INTERMEDIATE : $(CIRCOBJS) $(CIRCEQOBJS) $(IEOBJS) $(GKGOBJS) $(INTOBJS) $(LBOBJS) $(NROBJS) $(RRGWOBJS) $(SNTOBJS) $(SWSHOBJS)

#############################################################################
#
# All .cc files used to generate libraries.
#
#############################################################################

CIRCCC = CID.cc CKG.cc CKR.cc CTD.cc

CIRCEQCC = CEID.cc CEKG.cc CEKR.cc CETD.cc

IECC = IEKG.cc IEKR.cc IETD.cc

GKGCC = GKG.cc

INTCC = Integral.cc

LBCC = LB.cc

NRCC = 	NRElle.cc NREllf.cc NREllpi.cc NRFactrl.cc NRGammln.cc NRIndexx.cc \
        NRLUbksb.cc NRLUdcmp.cc NRPythag.cc NRUtil.cc \
        NRpzextr.cc NRrc.cc NRrd.cc NRrf.cc NRrj.cc \
        NRrzextr.cc NRspline.cc NRsplint.cc NRtqli.cc NRtred2.cc 

RRGWCC = RRGW.cc

SNTCC = SNTBSstep.cc SNTMmid.cc SNTOdeint.cc SNT.cc SNTSasNak.cc SNTTeuk.cc SNTInterp.cc

SWSHCC = SWSHCGUtil.cc SWSHSpherical.cc SWSHSpheroid.cc

#############################################################################
#
# Begin setting up dependencies.
#
#############################################################################

all : Circ Circ_Seq Circ_Traj Circ_Wave Circ_Eq Circ_Eq_Seq Circ_Eq_Traj \
Circ_Eq_Wave Circ_Eq_Rec GW Wrapper IE_test SchwarzRotate GenericSnapshot \
FastGenericSnapshot IEI

#############################################################################
#
# Dependencies for executables.
#
#############################################################################

Circ : Circ.cc -lCirc -lSWSH -lIntegral -lGKG -lSNT -lRRGW -lNR Globals.h CKG.h CKR.h CTD.h NRUtil.h RRGW.h SNT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ.cc -o $(BIN)/Circ -I$(INCL) -L$(LIB) -lCirc -lSWSH -lIntegral -lGKG -lSNT -lRRGW -lNR $(SYSLIBS)

Circ_Seq : Circ_Seq.cc -lCirc -lCircEq -lLB -lSWSH -lIntegral -lGKG -lSNT -lRRGW -lNR Globals.h CEKG.h CKG.h CEKR.h CKR.h CETD.h CTD.h NRUtil.h RRGW.h SNT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Seq.cc -o $(BIN)/Circ_Seq -I$(INCL) -L$(LIB) -lCirc -lCircEq -lLB -lSWSH -lIntegral -lGKG -lSNT -lRRGW -lNR $(SYSLIBS)

Circ_Traj : Circ_Traj.cc -lCirc -lLB -lGKG -lSWSH -lRRGW -lNR Globals.h CID.h NRUtil.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Traj.cc -o $(BIN)/Circ_Traj -I$(INCL) -L$(LIB) -lCirc -lLB -lGKG -lSWSH -lRRGW -lNR $(SYSLIBS)

Circ_Wave : Circ_Wave.cc -lCirc -lLB -lGKG -lSWSH -lRRGW -lNR Globals.h CID.h NRUtil.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Wave.cc -o $(BIN)/Circ_Wave -I$(INCL) -L$(LIB) -lCirc -lLB -lGKG -lSWSH -lRRGW -lNR $(SYSLIBS)

Circ_Eq : Circ_Eq.cc -lCircEq -lSWSH -lGKG -lSNT -lRRGW -lNR Globals.h CEKG.h CEKR.h NRUtil.h RRGW.h SNT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq.cc -o $(BIN)/Circ_Eq -I$(INCL) -L$(LIB) -lCircEq -lSWSH -lGKG -lSNT -lRRGW -lNR $(SYSLIBS)

Circ_Eq_Seq : Circ_Eq_Seq.cc -lCircEq -lSWSH -lGKG -lSNT -lRRGW -lNR Globals.h CEKG.h CEKR.h CETD.h NRUtil.h RRGW.h SNT.h SWSH.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Seq.cc -o $(BIN)/Circ_Eq_Seq -I$(INCL) -L$(LIB) -lCircEq -lSWSH -lGKG -lSNT -lRRGW -lNR $(SYSLIBS)

Circ_Eq_Traj : Circ_Eq_Traj.cc -lCircEq -lSWSH -lRRGW -lNR Globals.h CEID.h NRUtil.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Traj.cc -o $(BIN)/Circ_Eq_Traj -I$(INCL) -L$(LIB) -lCircEq -lSWSH -lRRGW -lNR $(SYSLIBS)

Circ_Eq_Wave : Circ_Eq_Wave.cc -lCircEq -lSWSH -lRRGW -lNR Globals.h NRUtil.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Wave.cc -o $(BIN)/Circ_Eq_Wave -I$(INCL) -L$(LIB) -lCircEq -lSWSH -lRRGW -lNR $(SYSLIBS)

Circ_Eq_Rec : Circ_Eq_Rec.cc -lCircEq -lSWSH -lRRGW -lNR GKG.h Globals.h NRUtil.h SWSH.h RRGW.h
	$(CC) $(CFLAGS) $(EXECSRC)/Circ_Eq_Rec.cc -o $(BIN)/Circ_Eq_Rec -I$(INCL) -L$(LIB) -lCircEq -lSWSH -lRRGW -lNR $(SYSLIBS)

GW : GW.cc -lSWSH -lRRGW -lNR Globals.h RRGW.h SWSH.h NRUtil.h
	$(CC) $(CFLAGS) $(EXECSRC)/GW.cc -o $(BIN)/GW -I$(INCL) -L$(LIB) -lSWSH -lRRGW -lNR $(SYSLIBS)

Wrapper : Wrapper.cc -lCirc -lIE -lSNT -lGKG -lIntegral -lSWSH -lNR Globals.h IEKG.h IEKR.h NRUtil.h NRIEKG.h 
	$(CC) $(CFLAGS) $(EXECSRC)/Wrapper.cc -o $(BIN)/Wrapper -I$(INCL) -L$(LIB) -lCirc -lIE -lSNT -lGKG -lIntegral -lSWSH -lNR $(SYSLIBS)

IE_test : IE_test.cc -lIE -lSNT -lGKG -lSWSH -lNR IEKG.h IEKR.h Globals.h SNT.h SWSH.h NRUtil.h NRIEKG.h  
	$(CC) $(CFLAGS) $(EXECSRC)/IE_test.cc -o $(BIN)/IE_test -I$(INCL) -L$(LIB) -lIE -lSNT -lGKG -lSWSH -lNR $(SYSLIBS)

SchwarzRotate : SchwarzRotate.cc -lIE -lSNT -lGKG -lSWSH -lNR IEKG.h IEKR.h Globals.h SNT.h SWSH.h NRUtil.h NRIEKG.h
	$(CC) $(CFLAGS) $(EXECSRC)/SchwarzRotate.cc -o $(BIN)/SchwarzRotate -I$(INCL) -L$(LIB) -lCirc -lIE -lSNT -lGKG -lIntegral -lSWSH -lNR $(SYSLIBS)

GenericSnapshot: GenericSnapshot.cc -lIE -lSNT -lGKG -lSWSH -lNR IEKG.h IEKR.h Globals.h SNT.h SWSH.h NRUtil.h NRIEKG.h
	$(CC) $(CFLAGS) $(EXECSRC)/GenericSnapshot.cc -o $(BIN)/GenericSnapshot -I$(INCL) -L$(LIB) -lIE -lSNT -lGKG -lSWSH -lNR $(SYSLIBS)

FastGenericSnapshot: FastGenericSnapshot.cc -lIE -lSNT -lGKG -lSWSH -lNR IEKG.h IEKR.h Globals.h SNT.h SWSH.h NRUtil.h NRIEKG.h
	$(CC) $(CFLAGS) $(EXECSRC)/FastGenericSnapshot.cc -o $(BIN)/FastGenericSnapshot -I$(INCL) -L$(LIB) -lIE -lSNT -lGKG -lSWSH -lNR $(SYSLIBS)

IEI: IEI.cc -lIE -lSNT -lGKG -lSWSH -lNR IEKG.h IEKR.h Globals.h SNT.h SWSH.h NRUtil.h NRIEKG.h
	$(CC) $(CFLAGS) $(EXECSRC)/IEI.cc -o $(BIN)/IEI -I$(INCL) -L$(LIB) -lIE -lSNT -lGKG -lSWSH -lNR $(SYSLIBS)


#############################################################################
#
# Dependencies for libraries.
#
#############################################################################

-lCirc : $(CIRCCC) $(CIRCOBJS) Globals.h CID.h CKG.h CKR.h CTD.h Integral.h NRCKG.h NRUtil.h SWSH.h SNT.h
	$(AR) $(LIB)/libCirc.a $(CIRCOBJS); \
	ranlib $(LIB)/libcirc.a;

-lCircEq : $(CIRCEQCC) $(CIRCEQOBJS) Globals.h CEID.h CEKG.h CEKR.h CETD.h NRUtil.h SWSH.h SNT.h
	$(AR) $(LIB)/libCircEq.a $(CIRCEQOBJS); \
	ranlib $(LIB)/libCircEq.a;

-lIE : $(IECC) $(IEOBJS) Globals.h IEKG.h IEKR.h IETD.h NRUtil.h NRIEKG.h SWSH.h SNT.h
	$(AR) $(LIB)/libIE.a $(IEOBJS); \
	ranlib $(LIB)/libIE.a;

-lGKG : $(GKGCC) $(GKGOBJS) Globals.h GKG.h NRUtil.h
	$(AR) $(LIB)/libGKG.a $(GKGOBJS); \
	ranlib $(LIB)/libGKG.a;

-lIntegral : $(INTCC) $(INTOBJS) Globals.h Integral.h NRUtil.h
	$(AR) $(LIB)/libIntegral.a $(INTOBJS); \
	ranlib $(LIB)/libIntegral.a;

-lLB : $(LBCC) $(LBOBJS) Globals.h GKG.h LB.h NRLB.h NRUtil.h
	$(AR) $(LIB)/libLB.a $(LBOBJS); \
	ranlib $(LIB)/libLB.a;

-lNR : $(NRCC) $(NROBJS) Globals.h NRUtil.h
	$(AR) $(LIB)/libNR.a $(NROBJS); \
	ranlib $(LIB)/libNR.a;

-lRRGW : $(RRGWCC) $(RRGWOBJS) Globals.h RRGW.h SWSH.h
	$(AR) $(LIB)/libRRGW.a $(RRGWOBJS); \
	ranlib $(LIB)/libRRGW.a;

-lSNT : $(SNTCC) $(SNTOBJS) Globals.h SNT.h NRSNT.h NRUtil.h
	$(AR) $(LIB)/libSNT.a $(SNTOBJS); \
	ranlib $(LIB)/libSNT.a;

-lSWSH : $(SWSHCC) $(SWSHOBJS) Globals.h SWSH.h NRSWSH.h
	$(AR) $(LIB)/libSWSH.a $(SWSHOBJS); \
	ranlib $(LIB)/libSWSH.a;

#############################################################################
#
# Rule for implicitly generating intermediate .o from .cc (needed to make
# libaries).
#
#############################################################################

%.o : %.cc
	$(CC) $(CFLAGS) -I$(INCL) -c $< -o $@

#############################################################################
#
# make clean
#
#############################################################################

clean : dummy
	$(RM) $(BIN)/*
	$(RM) $(LIB)/*
