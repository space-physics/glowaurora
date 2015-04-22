# USAGE: 
# make -s FC=<compiler>
# EXAMPLE:
# make -s FC=gfortran 
#

# use gfortran by default
ifeq ($(strip $(fc)),)
FC=gfortran
endif

FFLAGS = -O3 -I $(INCLUDE) -L $(LIBDIR)
#FFLAGS = -g $(DBGFLAGS) 

DBGFLAGS = -debug full -traceback
DBGFLAGS += -check bounds -check format -check output_conversion -check pointers -check uninit
DBGFLAGS += -fpe-all=0 # this traps all floating point exceptions

.SUFFIXES: .o .F .F90 .f

%.o: %.F90
	$(FC) $(FFLAGS) -c  -o $@ $<
%.o: %.F
	$(FC) $(FFLAGS) $(REAL8) -c  -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c  -o $@ $<
#
# Sources (in order of dependency):
#
SOURCES = aurexample.f ephoto.f egrid.f etrans.f exsect.f fieldm.f gchem.f geomag.f glow.f iri90.f maxt.f nrlmsise00.f qback.f rcolum.f rout.f snoem.f snoemint.f solzen.f ssflux.f vquart.f

OBJS := $(addsuffix .o, $(basename $(SOURCES)))
EXEC = auroraexample

$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

INCLUDE = /usr/include
LIBDIR = /usr/lib

clean:
	rm -f *.o *.mod $(EXEC)

