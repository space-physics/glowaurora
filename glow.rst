==============================================
Preliminary Documentation for the GLOW Model
==============================================

:Date: 4/05, 3/15
:Version: 0.973
:License: Open Source Academic Research License Agreement in file glowlicense.txt
:Author: Stan Solomon
:Phone: 303-497-2179
:Email: stans@ucar.edu
:Institution: HAO/NCAR, Boulder, CO 80307-3000
:Note: v0.973 is an minor incremental release from v0.97

**see README file for release notes**

.. contents::

Introduction
============

The GLOW family of subroutines calculate ionization
and excitation rates, energetic electron production and transport, excited
species densities, and airglow emission rates for the terrestrial thermosphere.
Use is governed by the Open Source Academic Research License Agreement
contained in the file ``glowlicense.txt``.  Standards and practices are
specified by the "rules of the road" adopted by the CEDAR and TIMED programs.
Other than the example driver programs, the header file glow.h, and
possibly the master subroutine glow.f, it is inadvisable to modify the code.
If any changes are required, please suggest them to Stan Solomon.

Since this is a developmental code, there are bound to be problems which
users will uncover.  Please let me know about them.  I will also attempt to
accommodate any reasonable suggestions for enhancements or changes.  By
following this procedure, users will benefit by staying compatible with future
developments, and from a more systematic elimination of any programming errors.
Also, please upgrade to the latest version.  Use of obsolete versions may lead
to incorrect results.

Program Description
===================

Programs DAYEXAMPLE and AUREXAMPLE are provided to give some guidance in
how to use the subroutine package for a daytime and auroral run, respectively.
Program HEXEXAMPLE shows how to do a high-energy electron flux calculation.
These programs use MSIS-2K and IRI-90 to specify the neutral atmosphere and
initial electron density profile.  The programs should be modified by the
user to suit particular purposes.  It is **not** necessary to use MSIS and IRI;
any model or measurement that specifies neutral densities, temperatures, and
high-altitude electron densities can be employed.  The output statements at the
end are just examples of a small subset of the calculated parameters.
Note that the user must change the directory specified in the call to IRI90
to correspond to the location of the ccir*.asc and ursi*.asc data files.

Subroutine GLOW is the master routine of the /glow package.  It receives
input parameters from the calling program in common block /CGLOW/, calls
the other subroutines, and returns results in /CGLOW/.  Header file glow.h
supplies array sizes for the altitude and electron energy grids.  It is
not necessary to use subroutine GLOW to call the others - those particularly
bold and well-informed users who desire to use some subset of the package
may attempt to do so, using subroutine GLOW as guidance.

The header file glow.h is used to specify the number of altitude levels,
energy bins, and wavelength bins.  It may be edited to select numbers
appropriate for the altitude and energy range under study before compiling.
The number of altitudes must correspond to the altitude array supplied by
the driver program.

The data files ``ephoto_x*.dat`` and ``ssflux_*.dat`` must exist on the **current
working directory**.  The number of wavelenth bins in these files must equal
the parameter LMAX specified in glow.h, and the wavelength ranges must be the
same in all input files.

These routines are written in standard Fortran-77, and should be
compatible with most compilers, including Fortran-90 compilers.  However,
their performace has only been tested using Portland Group Fortran-77 (pgf77)
and Intel Fortran-90 (ifort) compilers running under Linux Centos 6.3.  I
generally compile using third-level optimization (-O3).  I would be interested
in hearing of experiences with other compilers, but make no claims or promises.

Subroutine description
======================

Subroutines called by the example programs include:

==========  =============
Subroutine  Description
==========  =============
GTD7        NRL Mass Spectrometer Incoherent Scatter (MSIS-2K) model
TSELEC      Set-up entry for MSIS-2K
IRI90       International Reference Ionosphere-1990 (IRI-90) model (Belitza, 1990)
GEOMAG      translates geographic to geomagnetic coordinates and vice versa
MAXT        generates a Maxwellian electron spectrum
SNOEMINT    interpolates nitric oxide profile from SNOEM
SNOEM       special version of NOEM (Marsh et al., 2004)
==========  =============

Subroutines called by GLOW are:

==========  =============
Subroutine  Description
==========  =============
  EGRID     sets up electron energy grid
  FIELDM    calculates magnetic dip angle
  SOLZEN    calculates solar zenith angle
  SSFLUX    scales solar flux for activity level
  RCOLUM    calculates slant column density of major species
  EPHOTO    calculates photoionization and photoelectron production
  QBACK     estimates background ionization
  ETRANS    computes electron transport, ionization, excitation using Nagy & Banks 2-stream method calls EXSECT for cross-sections, first call only
  GCHEM     finds electron/ion/metastable densities, airglow emissions uses VQUART to solve electron density equation
==========  =============

Definitions
===========
(taken from the introductory comment to subroutine GLOW)

Supplied to subroutine in labeled common /CGLOW/:

==========  =============
Variable    Description
==========  =============
IDATE       Date, in form yyddd
UT          Universal Time; seconds
GLAT        Geographic latitude; degrees
GLONG       Geographic longitude; degrees
ISCALE      Solar flux scaling switch, see subroutine SSFLUX
JLOCAL      =0 for electron transport calculation, =1 for local calc only
KCHEM       Ion/electron chemistry switch, see subroutine GCHEM
F107        Solar 10.7 cm flux for day being modeled, 1.E-22 W m-2 Hz-1
F107A       Solar 10.7 cm flux 81-day centered average
HLYBR       H Ly-b (1026A) enhancement ratio (optional)
FEXVIR      Fe XVI (335A) enhancement ratio (optional)
HLYA        H Ly-a flux; photons cm-2 s-1 (optional)
HEIEW       He I 10830 equivalent width; milliAngstroms (obsolete)
XUVFAC      Factor by which to multiply solar flux 18-250A or 18-50A (optional)
ZZ          altitude array; cm
ZO          O number density at each altitude; cm-3
ZN2         N2  "      "      "   "     "       "
ZO2         O2         "
ZNO         NO         "
ZNS         N(4S)      "
ZND         N(2D)      "
ZRHO        mass density at each altitude; gm cm-3
ZE          electron density at each alt; cm-3
ZTN         neutral temperature at each alt; K
ZTI         ion temperature at each alt; K
ZTE         electron temp at each alt; K
PHITOP      energetic electron flux into top of atmosphere; cm-2 s-1 eV-1
EFLUX       (obsolete)
==========  =============

Calculated by subroutine:

==========  =============
Subroutine  Calculates
==========  =============
 SZA        solar zenith angle; radians
 DIP        magnetic field dip angle; radians
 EFRAC      energy conservation check from ETRANS, (out-in)/in
 IERR       error code returned from ETRANS:
            0=normal, 1=local problem, 2=transport problem
 ZMAJ       major species density array, O, O2, N2; cm-3
 ZCOL       major species slant column density array, O, O2, N2; cm-2
 WAVE1      longwave edge of solar flux wavelength range; A
 WAVE2      shortwave edge of solar flux wavelength range; A
 SFLUX      scaled solar flux in each wavelength range; photons cm-2 s-1
 ENER       electron energy grid; eV
 DEL        width of each bin in electron energy grid; eV
 PESPEC     photoelectron production rate at energy, altitude; cm-3 s-1
 SESPEC     secondary electron production rate (from EAURI); cm-3 s-1
 PHOTOI     photoionization rates for state, species, altitude; cm-3 s-1
            O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
            O2+ states: X, a+A, b, dissoc.
            N2+ states: X, A, B, C, F, dissoc.
 PHOTOD     photodissoc. & exc. rates for state, species, alt.; cm-3 s-1
            (1,2,J) = O2 -> O(3P) + O(1D))
            (2,2,J) = O2 -> O(3P) + O(1S)
            (1,3,J) = N2 -> N + N
 PHONO      photoionization/dissociation/excitation rates for NO, cm-3 s-1
            (1,J) = NO+ from H Ly-a ionization
 QTI        (obsolete)
 AURI       (obsolete)
 PIA        (obsolete)
 SION       electron impact ioniz. rates calculated by ETRANS; cm-3 s-1
 UFLX       upward hemispherical electron flux; cm-2 s-1 eV-1
 DFLX       downward hemispherical electron flux; cm-2 s-1 eV-1
 AGLW       Electron impact exc. rates; state, species, alt.; cm-3 s-1
            O states: 1D, 1S, 5S, 3S, 3p5P, 3p3P, 3d3D, 3s'3D
            O2 states: a, b, (A+A'+c), B(SRC), 9.9eV, Ryds., vib.
            N2 states: (A+B+W), B', C, (a+a'+w), 1Pu, b', Ryds., vib.
 EHEAT      ambient electron heating rate, eV cm-3 s-1
 TEZ        total energetic electron energy deposition, eV cm-3 s-1
 ECALC      electron density, calculated below 200 km, cm-3
 ZXDEN      array of excited and and/or ionized state densities: O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P), N(2D), O(1S), O(1D), 8 spares, at each altitude; cm-3
 ZETA       array of volume emission rates: 3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A, 7774A, 8446A, 3726A, 9 spares; cm-3 s-1
 ZCETA      array of contributions to each v.e.r at each alt; cm-3 s-1
 VCB        array of vertical column brightnesses (as above); Rayleighs
==========  =============

==========  =============
Array       Dimension (length):
==========  =============
JMAX        number of altitude levels
NBINS       number of energetic electron energy bins
LMAX        number of wavelength intervals for solar flux
NMAJ        number of major species
NEX         number of ionized/excited species
NW          number of airglow emission wavelengths
NC          number of component production terms for each emission
NST         number of states produced by photoionization/dissociation
NEI         number of states produced by electron impact
NF          (obsolete)
==========  =============

Additional quantites (have not verified definitions M.H. Aug 2015)

==========  =============
Variable    Description
==========  =============
TeFlux      Total electron flux ?

==========  =============

Notes:
======

MSIS-2000 and IRI-90
--------------------
Versions of MSIS-2K and IRI-90 are provided for the convenience of users
who do not have their own copies.  Attribution to the appropriate sources
(Hedin, 1991; Picone, 2002; Belitza, 1990) should be made.  This is not the
standard version of IRI, as I have modified it to make it work on various
systems, particularly volatile memory systems.  However, I cannot guarantee
that the results obtained from them are correct.  Caution: IRI occasionally
writes mysterious messages on unit 12.  The MSIS-2K subroutine (GTD7) may
be obtained by from http://download.hao.ucar.edu/pub/stans/msis.
The IRI90 subroutine in file iri90.f and its data files ccir*.asc and ursi*.asc
may be obtained by from http://download.hao.ucar.edu/pub/stans/iri.
Note that the user must change the directory specified in the call to IRI90
to correspond to the location of the ccir*.asc and ursi*.asc data files.

glow.h
------
The header file ``glow.h`` is "included" in the routines that require the
altitude, electron energy, and wavelength grid size parameters.  Execution
time increases rapidly with size of the energy grid.  The default grid extends
to 50 keV and is appropriate for most photoelectron and auroral calculations;
the electron energy grid can extend (in principle) to 1 GeV.  The altitude
grid may be altered to suit the user's needs, but with caution.  The two
biggest pitfalls are not providing enough resolution and not providing a deep
enough atmosphere.  There should be no significant flux of electrons out of
the bottom of the altitude grid.  As for resolution, a rule of thumb is about
four points per scale height.

Sanity checks
-------------
The ETRANS error code IERR and energy conservation ratio EFRAC
should be checked for normal return.  IERR should equal zero; if it doesn't
it means that the total inelastic cross section is near zero somewhere, which
is usually caused by a near zero ambient electron density.  EFRAC should be
less than ~0.03 for photoelectron calculations, and less than ~0.1 for aurora.
It can get up to the 0.1-0.2 range in the twilight, which is not good, but at
present unavoidable.

Ne caveats
----------
Electron density calculations can be made by GCHEM below 200 km but not
above where transport/diffusion effects become important.  Therefore, an
electron density profile (such as from the NCAR TIE-GCM or from IRI) must be
provided above 200 km.  An initial non-zero electron density profile must be
provided at all altitudes in array ZE because otherwise ETRANS will produce
an error.  Calculated electron densities are returned in array ECALC, with
values from ZE included where calculations are not made.  The switch KCHEM
determines what ion/electron density calculations are made (see subroutine
GCHEM).  With the exception of KCHEM=0 (no calculations), none are foolproof.
One irritation that may be noted when KCHEM=4 (calculate Ne below 200 but use
provided Ne above 200) is a discontinuity in Ne(z) at 200 km.  GCHEM attempts
to deal with this by interpolating from the 200 km level to three grid
points above it.  Another caveat with KCHEM=4 or KCHEM=3 (Ne provided,
ions calculated) is that if there is an incompatibility between the specifed
Ne(z) and the calculated ionization rates, negative values for some ions may
result.  KCHEM=2 (electrons and major ions provided, minor ions calculated) and
KCHEM=1 (electrons and all ions except O+(2D,2P) provided) should be fairly
reliable.

XUV
---
For daytime calculations, the parameter XUVFAC is provided to deal with the
uncertainty concerning the solar spectrum from 18-250 A.  When the
Hinteregger et al. [1981] model is employed, a reasonable value for XUVFAC is
2.0, as suggested by Richards et al. [1984; 1994], however a value as high as
4.0 might be realistic, as discussed by Solomon et al. [2001].  When the
EUVAC model [Richards et al., 1994] is employed, XUVFAC is only applied
to the region still obtained from the Hinteregger model, 18-50 A, since
longward of that point EUVAC has already increased the solar fluxes relative
to the Hinteregger spectrum (by factors of 2-3).  SSFLUX also now provides
the ability for the user to specify a solar spectrum from other models or
measurements, in which case XUVFAC is ignored.  The number of bins in the
solar spectrum input file must be equal to LMAX (specified in glow.h) and
the wavelength ranges must correspond to those in ephoto_x*.dat.

Low energy precipitation
------------------------
The upper altitude boundary of the electron transport calculation by ETRANS is
specified by the PHITOP array, which may contain a flux of auroral electrons,
conjugate photoelectrons, or both.  In the program AUREXAMPLE, an initial
electron density profile is obtained from IRI for the first call to GLOW,
then it is replaced by the calculated profile below 200 km (and a constant
value above 200 km), and GLOW is called again.  This isn't really necessary
but it gives an improved estimate of the low-energy electron flux (which
depends on the ambient electron density), since IRI is not valid in the auroral
regions.  For high-energy calculations this second call may safely be skipped.
In any case, it is safest to put a floor on the electron density profile, e.g., ZE(J) > 100.

NO Density
----------
The example programs contain an estimate of NO density from the NOEM
empirical model (Marsh et al., 2004), which is based on measurements by the
SNOE satellite.  This can be important for the NO+/O2+ ratio in the lower
ionosphere and hence has a small effect on E-region electron density, but does
not otherwise significantly affect the model.

Electron Impact Cross-section
-----------------------------
Electron impact cross sections employed by the model can be obtained from
common block CXSECT if necessary; just include a copy of this common block
(from EXSECT) in the calling program.

Bibliography
============
.. [1] Nagy, A. F., and P. M. Banks, Photoelectron fluxes in the ionosphere, J. Geophys. Res., 75, 6260, 1970.
.. [2] Solomon, S. C., P. B. Hays, and V. J. Abreu, The auroral 6300A emission: Observations and modeling, J. Geophys. Res., 93, 9867, 1988.
.. [3] Solomon, S. C., and V. J. Abreu, The 630 nm dayglow, J. Geophys. Res., 94, 6817, 1989.
.. [4] Solomon, S. C., Auroral particle transport using Monte Carlo and hybrid  methods, J. Geophys. Res., 106, 107, 2001.
.. [5] Solomon, S. C., S. M. Bailey, and T. N. Woods, Effect of solar soft X-rays on the lower atmosphere, Geophys. Res. Lett., 28, 2149, 2001.
.. [6] Bailey, S. M., C. A. Barth, and S. C. Solomon, A model of nitric oxide in the lower thermosphere, J. Geophys. Res., 107, 1205, 2002.
