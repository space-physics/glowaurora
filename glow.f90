! Subroutine GLOW
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Version 0.973
!
! Stan Solomon, 1988, 1989, 1990, 1991, 1992, 1994, 2000, 2002, 2005, 2015
!
! Subroutine GLOW is the master routine of the /glow package.  It
! receives input parameters from the calling program in common block
! /CGLOW/, calls the other subroutines, and returns results in /CGLOW/.
! Header file glow.h supplies array sizes for the altitude, electron,
! and solar spectrum energy grids.
!
! Subroutines called by GLOW are:
!   EGRID   sets up electron energy grid
!   FIELDM  calculates magnetic dip angle
!   SOLZEN  calculates solar zenith angle
!   SSFLUX  scales solar flux for activity level
!   RCOLUM  calculates slant column density of major species
!   EPHOTO  calculates photoionization and photoelectron production
!   QBACK   estimates background ionization
!   ETRANS  computes electron transport, ionization, excitation
!             calls EXSECT for cross-sections, first call only
!   GCHEM   finds electron/ion/metastable densities, airglow emissions
!             uses VQUART to solve electron density equation
!
! Supplied to subroutine in labeled common /CGLOW/:
! IDATE   Date, in form yyddd
! UT      Universal Time; seconds
! GLAT    Geographic latitude; degrees
! GLONG   Geographic longitude; degrees
! ISCALE  Solar flux scaling switch, see subroutine SSFLUX
! JLOCAL  =0 for electron transport calculation, =1 for local calc only
! KCHEM   Ion/electron chemistry switch, see subroutine GCHEM
! F107    Solar 10.7 cm flux for day being modeled, 1.E-22 W m-2 Hz-1
! F107A   Solar 10.7 cm flux 81-day centered average
! HLYBR   H Ly-b (1026A) enhancement ratio
! FEXVIR  Fe XVI (335A) enhancement ratio
! HLYA    H Ly-a flux; photons cm-2 s-1
! HEIEW   He I 10830 equivalent width; milliAngstroms (obsolete)
! XUVFAC  Factor by which to multiply to solar flux 16-250 A or 16-50 A.
! ZZ      altitude array; cm
! ZO      O number density at each altitude; cm-3
! ZN2     N2  "      "      "   "     "       "
! ZO2     O2         "
! ZNO     NO         "
! ZNS     N(4S)      "
! ZND     N(2D)      "
! ZRHO    mass density at each altitude; gm cm-3
! ZE      electron density at each alt; cm-3
! ZTN     neutral temperature at each alt; K
! ZTI     ion temperature at each alt; K
! ZTE     electron temp at each alt; K
! PHITOP  energetic electron flux into top of atmosphere; cm-2 s-1 eV-1
! EFLUX   obsolete
! EZERO   obsolete
!
! Calculated by subroutine:
! SZA     solar zenith angle; radians
! DIP     magnetic field dip angle; radians
! EFRAC   energy conservation check from ETRANS, (out-in)/in
! IERR    error code returned from ETRANS:
!           0=normal, 1=local problem, 2=transport problem
! ZMAJ    major species density array, O, O2, N2; cm-3
! ZCOL    major species slant column density array, O, O2, N2; cm-2
! WAVE1   longwave edge of solar flux wavelength range; A
! WAVE2   shortwave edge of solar flux wavelength range; A
! SFLUX   scaled solar flux in each wavelength range; photons cm-2 s-1
! ENER    electron energy grid; eV
! DEL     width of each bin in electron energy grid; eV
! PESPEC  photoelectron production rate at energy, altitude; cm-3 s-1
! SESPEC  background electron production rate (obsolete); cm-3 s-1
! PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
!           O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
!           O2+ states: X, a+A, b, dissoc.
!           N2+ states: X, A, B, C, F, dissoc.
! PHOTOD  photodissoc. & exc. rates for state, species, alt.; cm-3 s-1
!           (1,2,J) = O2 -> O(3P) + O(1D))
!           (2,2,J) = O2 -> O(3P) + O(1S)
!           (1,3,J) = N2 -> N + N
! PHONO   photoionization/dissociation/excitation rates for NO, cm-3 s-1
!         (1,J) = NO+ from H Ly-a ionization
! QTI     obsolete
! AURI    obsolete
! PIA     obsolete
! SION    electron impact ioniz. rates calculated by ETRANS; cm-3 s-1
! UFLX    upward hemispherical electron flux; cm-2 s-1 eV-1
! DFLX    downward hemispherical electron flux; cm-2 s-1 eV-1
! AGLW    Electron impact exc. rates; state, species, alt.; cm-3 s-1
!           O states: 1D, 1S, 5S, 3S, 3p5P, 3p3P, 3d3D, 3s'3D
!           O2 states: a, b, (A+A'+c), B(SRC), 9.9eV, Ryds., vib.
!           N2 states: (A+B+W), B', C, (a+a'+w), 1Pu, b', Ryds., vib.
! EHEAT   ambient electron heating rate, eV cm-3 s-1
! TEZ     total energetic electron energy deposition, eV cm-3 s-1
! ECALC   electron density, calculated below 200 km, cm-3
! ZXDEN   array of excited and and/or ionized state densities:
!           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
!           N(2D), O(1S), O(1D), 8 spares, at each altitude; cm-3
! ZETA    array of volume emission rates:
!           3371A, 4278A, 5200A, 5577A, 6300A, 7320A, 10400A, 3466A,
!           7774A, 8446A, 3726A, 9 spares; cm-3 s-1
! ZCETA   array of contributions to each v.e.r at each alt; cm-3 s-1
! VCB     array of vertical column brightnesses (as above); Rayleighs
!
! Array dimensions:
! JMAX    number of altitude levels
! NBINS   number of energetic electron energy bins
! LMAX    number of wavelength intervals for solar flux
! NMAJ    number of major species
! NEX     number of ionized/excited species
! NW      number of airglow emission wavelengths
! NC      number of component production terms for each emission
! NST     number of states produced by photoionization/dissociation
! NEI     number of states produced by electron impact
! NF      obsolete
!
!

module glowmod
    use ccglow
    private
    public :: glow

contains
      SUBROUTINE GLOW
      use rcolummod
      use szacalc
      use qbackmod
!
      PARAMETER (NEX=20)
      PARAMETER (NW=20)
      PARAMETER (NC=10)
      PARAMETER (NF=4)

      Real(dp) :: ENER, DEL
!
      COMMON /CGLOW/ &
         IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM, &
         F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC, &
         ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX), &
         ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX), &
         ZTN(JMAX), ZTI(JMAX), ZTE(JMAX), &
         PHITOP(NBINS), EFLUX(NF), EZERO(NF), &
         SZA, DIP, EFRAC, IERR, &
         ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX), &
         WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX), &
         ENER(NBINS), DEL(NBINS), &
         PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX), &
         PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX), &
         QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX), &
         UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX), &
         EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX), &
         ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW)
!
      DIMENSION ZVCD(NMAJ,JMAX)

      REAL(sp) TEFLUX
!
!      DATA IFIRST/1/
!
!
! First call only: set up energy grid:
!
!      IF (IFIRST .EQ. 1) THEN
!        IFIRST = 0
!        CALL EGRID (ENER, DEL, NBINS)
!      ENDIF
!
!
! Find magnetic dip angle and solar zenith angle (radians):
!
      CALL FIELDM (GLAT, GLONG, 300., XF, YF, ZF, FF, DIP, DEC, SDIP)
      DIP = ABS(DIP) * PI/180.
!
      CALL SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
      SZA = SZA * PI/180.
!
!
! Scale solar flux:
!
    CALL SSFLUX (ISCALE, F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC, WAVE1, WAVE2, SFLUX)
!
!
! Pack major species density array:
!
      DO J=1,JMAX
        ZMAJ(1,J) = ZO(J)
        ZMAJ(2,J) = ZO2(J)
        ZMAJ(3,J) = ZN2(J)
      End Do
!
!
! Calculate slant path column densities of major species in the
! direction of the sun:
!
        CALL RCOLUM (SZA, ZZ, ZMAJ, ZTN, ZCOL, ZVCD, NMAJ)
!
!
! Call subroutine EPHOTO to calculate the photoelectron production
! spectrum and photoionization rates as a function of altitude,
! unless all altitudes are dark, in which case zero arrays:
!
      IF (SZA .LT. 2.) THEN
        CALL EPHOTO
      ELSE
        PHOTOI = 0.
        PHOTOD = 0.
        PHONO = 0.
        PESPEC = 0.
      ENDIF
!
!
! Zero obsolete auroral primary and secondary arrays:
!
      PIA = 0.0
      SESPEC = 0.0
!
! Add background ionization to photoionization:
!
      CALL QBACK (ZMAJ, ZNO, ZVCD, PHOTOI, PHONO)
!
!
! Call subroutine ETRANS to calculate photoelectron and auroral
! electron transport and electron impact excitation rates, unless
! there are no energetic electrons, in which case zero arrays:
!
      TEFLUX = 0.
      DO N=1,NBINS !no vectorize
        TEFLUX = TEFLUX + PHITOP(N)
      End DO
!
      IF (TEFLUX.GT.0.001 .OR. SZA.LT.2.) THEN
        CALL ETRANS
      ELSE
        UFLX = 0.
        DFLX = 0.
        SION = 0.
        AGLW = 0.
        EHEAT = 0.
        TEZ = 0.
        EFRAC = 0.
        IERR = 0
      ENDIF
!
!
! Call subroutine GCHEM to calculate the densities of excited and
! ionized consituents, airglow emission rates, and vertical column
! brightnesses:
!
      CALL GCHEM
!
!
      END Subroutine GLOW
end module glowmod
