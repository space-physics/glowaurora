# USAGE: 
# make -s FC=<compiler>
# EXAMPLE:
# make -s FC=gfortran 
#

# use gfortran by default
ifeq ($(strip $(fc)),)
FC=gfortran
endif

FFLAGS = -O3 -I $(INCLUDE) -L $(LIBDIR) -fno-align-commons 
#-fdefault-real-8 SEG FAULT
# -ffpe-trap=precision -g
#FFLAGS = -g $(DBGFLAGS) 

DBGFLAGS = -debug full -traceback
DBGFLAGS += -check bounds -check format -check output_conversion -check pointers -check uninit
DBGFLAGS += -fpe-all=0 # this traps all floating point exceptions

.SUFFIXES: .o .F .F90 .f90 .f .mod

%.o: %.F90
	$(FC) $(FFLAGS) -c  -o $@ $<
%.o: %.f90
	$(FC) $(FFLAGS) -c  -o $@ $<
%.o: %.F
	$(FC) $(FFLAGS) $(REAL8) -c  -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c  -o $@ $<
#
# Sources (in order of dependency):
#
SOURCES = egrid.f ephoto.f maxt.f90 rcolum.f etrans.f exsect.f fieldm.f vquart.f90 gchem.f geomag.f solzen.f90 glow.f iri90.f nrlmsise00.f qback.f rout.f snoem.f90 snoemint.f ssflux.f aurexample.f 

OBJS := $(addsuffix .o, $(basename $(SOURCES)))
EXEC = auroraexample

$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

INCLUDE = /usr/include
LIBDIR = /usr/lib

clean:
	rm -f *.o *.mod $(EXEC)

