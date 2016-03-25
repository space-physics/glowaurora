# USAGE: 
# make -s FC=<compiler>
# EXAMPLE:
# make -s FC=gfortran 

AURORA=auroraexample

HEX=hexexample

PYMOD=glowfort

SP=./fortran/
#====== COMPILERS ==============================================================
# use gfortran by default -- Gfortran >= 5 required!
ifeq ($(strip $(fc)),)
FC=gfortran
endif

ifeq ($(FC),gfortran)
    #GFORTRAN >= 5
    FC_MAJOR=$(shell gfortran -dumpversion | cut -f1 -d.)
    ifeq ($(FC_MAJOR),4)
    $(warning "must have Gfortran >= 5")
    endif
    $(info $(FC_MAJOR))
    FFLAGS = -O3 -march=native -fno-align-commons
    #FFLAGS += -Wall -Wextra -fbounds-check
    #FFLAGS += -ffpe-trap=precision -g -ftraceback
    #FFLAGS += -fdefault-real-8 SEG FAULT
else ifeq ($(FC),ifort)
    #IFORT
    FFLAGS = -O3 -march=native -fno-align-commons
    #FFLAGS += -debug full -traceback  -C -fpe:0
endif

F2PY = f2py3
PYFLAGS = --quiet
#PYFLAGS += --f90exec=gfortran-5
#PYFLAGS =  --opt='-fno-align-commons'
#PYFLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1
#PYFLAGS +=--noarch --f90flags='${ARCHFLAGS}'
#===============================================================================

.SUFFIXES: .o .f90 .f .mod

%.o: %.f90
	$(FC) $(FFLAGS) -c  -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c  -o $@ $<
#
# Sources (in order of dependency):
#
FSRC = egrid.f ephoto.f maxt.f rcolum.f etrans.f exsect.f fieldm.f vquart.f gchem.f geomag.f solzen.f glow.f iri90.f nrlmsise00.f qback.f rout.f snoem.f snoemint.f ssflux.f

FAUR = aurexample.f
FHEX = hexexample.f

FSUB = aurora_sub.f

FLIBS = -latlas -llapack -lblas -lpthread


PSRC = $(addprefix $(SP),$(FSRC))
PAUR = $(addprefix $(SP),$(FAUR))
PHEX = $(addprefix $(SP),$(FHEX))
PSUB  = $(addprefix $(SP),$(FSUB))

FOBJ = $(PSRC:.f=.o)
AOBJ = $(PAUR:.f=.o)
HOBJ = $(PHEX:.f=.o)
#=====================================================

all: $(AURORA) $(HEX) $(PYMOD)

aurora: $(AURORA)

hex: $(HEX)

py: $(PYMOD)

$(AURORA): $(FOBJ) $(AOBJ)
	$(FC) $(FFLAGS) -o $@ $(FOBJ) $(AOBJ) $(FLIBS) $(LDFLAGS)

$(HEX): $(FOBJ) $(HOBJ)
	$(FC) $(FFLAGS) -o $@ $(FOBJ) $(HOBJ) $(FLIBS) $(LDFLAGS)

$(PYMOD): $(PSRC) $(PSUB)
	$(F2PY) $(PYFLAGS) -m $@ -c $(PSRC) $(PSUB) $(FLIBS) $(LDFLAGS)

clean:
	$(RM) $(FOBJ) $(MOBJ)

#printout a make variable by make print-myvariable
print-%  : ; @echo $* = $($*)

