
Z620=1
#-------------------------------------------------------------------------------
ifdef Z620
  CPP    = fpp
  FC     = ifort
  CC     = 
  NCDIR  = /usr/local/netcdf
  DFLAG  = -lnetcdf
  FFLAGS = -vec_report0
  SLIBS  = 
endif
ifdef SW
  CPP    = cpp
  FC     = mpif90
  CC     = mpicc
  NCDIR  = /public/software/mathlib/netcdf401
  DFLAG  = #-lnetcdf
  FFLAGS = -check all -i4 -r8  -g
  SLIBS  = #-lmpi
endif
ifdef SZ
  CPP    = cpp
  FC     = mpiswf90
  CC     = mpiswcc
  NCDIR  = /home/export/base/fiogrp/local_sz3
  DFLAG  = -DCP_SWF90 -L$(NCDIR)/lib  -lnetcdf
  #FFLAGS = -check_bounds -i4 -r8  -g -OPT:IEEE_arith=1 -I$(NCDIR)/include
  FFLAGS = -i4 -r8 -O2  -OPT:IEEE_arith=1 -I$(NCDIR)/include
  SLIBS  = -lsw
endif

#-------------------------------------------------------------------------------

INCLDIR    := -I. -I$(NCDIR)/include
SLIBS      := -L${NCDIR}/lib -lnetcdf $(SLIBS)

VPATH = ../
INC   =
EXE   = globalrun

SRC   = netcdf_mod.f90 time_mod.f90 masnum3d.f90

#-------------------------------------------------------------------------------

CPP1 = 
CPPE = 
# Function: -P : Inhibit generation of linemarkers.
#           -w : Avoid warning due to undefine __FILE__ and __LINE__.
#           -U__FILE__ -U__LINE__ : undefine __FILE__ and __LINE__ ,for use in future.

SEDP = 
# Function: Add char '~' befaore  '\' at end of line
#           Replace '//' to '/\' ,for '//' mean concat string in fortran,but mean comment in c
#    where: \ / is special char in sed command line,use \\,\/.
#           $ is special char in makefile,use $$ .
#           $ means end of line at command of sed.


SEDN = 
# Function: (1) Replace ";~" and any spaces into a newline
#           (2) Replace "~" into "&\n"

SEDE = 
#           (1) Delete characters before & for the case that there
#               is none valid characters before this &.
#           (2) Replace '/\' back to '//'

#-------------------------------------------------------------------------------
all : objs
	echo $0
	cd objs;pwd;make -f ../makefile   $(EXE);cp $(EXE) ..
#all : $(EXE)
objs :
	mkdir objs

OBJS = ${SRC:.f90=.o}

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

# Function: 1st step expand src; 2nd step compile src
.f90.o:
	$(FC) -c $(INCLDIR) $(FFLAGS) $<

$(EXE)  : $(OBJS)
	$(FC) -o $@  $(OBJS) $(SLIBS)

#-------------------------------------------------------------------------------

# Only make expanded src
#.F90.f90:
#	$(SEDP) $< | $(CPP1) | $(SEDN)  |$(CPPE) -D__FILE__=\"$*.f90\" | $(SEDE)>$*.f90
#	$(SEDP) $< | $(CPP1) | $(SEDE) >$*.f90

# clean expanded src
cleanesrc:
	rm $(esrc);	echo>cleanesrc

# clean expanded src first
expandsrc: cleanesrc $(esrc)
	rm cleanesrc

#-------------------------------------------------------------------------------

clean:
	cd objs;pwd;rm -f  $(EXE) *.o *.mod

cleano:
	rm -f  *.o *.mod

#-------------------------------------------------------------------------------

#$(OBJS) : $(INC)

#irrp_kernal_mod.o: irrp_smpi_mod.o irrp_split_mod.o
#irrp_package_mod.o: irrp_kernal_mod.o irrp_split_mod.o
#def_common_mod.o: time_mod.o netcdf_mod.o
#masnump1.o: def_common_mod.o

#-------------------------------------------------------------------------------

