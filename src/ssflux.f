! Subroutine SSFLUX
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Subroutine SSFLUX calculates the solar EUV and FUV flux in the range
! 0.5 to 1750 Angstroms for a specified level of solar activity.
!
! The calling routine supplies a scaling switch ISCALE, the daily 10.7
! cm flux F107, its 81-day centered average F107A, and, optionally,
! the H Lyman-Beta 1026A ratio to its solar minimum value HLYBR, the
! Fe XVI 335A ratio to its solar minimum value FEXVIR, the H Lyman-alpha
! 1216A flux HLYA, the He I 10830A equivalent width HEIEW, and an XUV
! enhancement factor XUVFAC.  Any optional calling parameters not used
! should be set to zero.
!
! XUVFAC is applied from 18-250 A for the Hinteregger model (ISCALE=0),
! from 18-50 A for the EUVAC model (ISCALE=1), and not at all for
! user-supplied data (ISCALE=2)
!
! The subroutine returns the longwave boundary WAVE1 and shortwave
! boundary WAVE2 of the wavelenth bins, and the solar flux in each bin
! SFLUX.  Bins are arranged in energy multiples of 2 from 0.5 to 8 A,
! aligned with k-shell boundaries at 23, 32, and 44 A from 18 to 44 nm,
! 10 A in width from 60 to 1050 A, and 50 A in width from 1050 to
! 1750 A with the exception of Lyman-alpha which has its own bin from
! 1210 to 1220 A.
!
! Methods used:
!   If ISCALE=0 the flux is scaled using parameterization methods based
! on F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio
! method (Hinteregger et al., GRL, 8, 1147, 1981) is used, based on the
! reference spectrum SC#21REFW at 1 nm bin resolution.  If the
! H Lyman-Beta (1026A) or Fe XVI (335A) enhancement ratios are provided
! (>0) as calling arguments, they are used to scale the EUV spectrum.
! Otherwise, enhancement ratios for H Ly-B and Fe XVI are calculated
! from F107 and F107A using Hinteregger's formula, employing
! coefficients which reduce to the reference values at F107=67.6,
! F107A=71.5.  The 'best fit' coefficients are not used as they produce
! some negative values at low solar activity, but remain in a
! 'commented out' data statement for reference.  The EUV spectrum is
! then scaled from these modeled ratios.  Scaling factors were
! calculated from contrast ratios in the SC#21REFW data file.
!   If ISCALE=1, the EUV flux (50-1050A) is scaled using the EUVAC model
! (Richards et al., JGR 99, 8981, 1994) re-binned onto ~1 nm intervals.
! The Hinteregger spectrum, scaled using the EUVAC algorithm, is used
! from 18 to 50A.
!   Neither of these models extends shortward of 18A, so from 1-18 A
! an amalgam of sources are used to derive an estimated flux, e.g.,
! DeJager, in Astronomical Observations from Space Vehicles, Steinberg,
! ed., 1964; Smith & Gottlieb, SSR 16, 771, 1974; Manson, in The Solar
! Output and its Variation, White, ed., 1977; Kreplin et al, ibid;
! Horan & Kreplin, Solar Physics 74, 265, 1981; Wagner, Adv. Space Res.
! 8, (7)67, 1988.
!    For FUV from 1050A-1750A, 50A interval bins from the Woods and
! Rottman [2002] reference spectrum and scale factors based on
! UARS SOLSTICE data are used.  The scaling method follows the
! Hinteregger or EUVAC algorithm, whichever is selected, so as to
! linearly scale the spectrum between the reference value and maximum
! value calculated with F10.7=F10.7A=200.  If a value for Lyman-alpha
! (HLYA>0) is provided by the calling program, it is subsituted into
! the spectrum.
!   If ISCALE=2, the solar flux (0-1750A) is read from a file named
! ssflux_user.dat in the current working directory.  The file must
! contain three columns:  WAVES, WAVEL, SFLUX (Angstroms and cm-2 s-1)
! in order of increasing wavelength.  The number of lines in the file
! must match the value of LMAX in glow.h.
!
! Modification history:
!   Stan Solomon, 12/88  Basic Hinteregger EUV, approx. SME FUV
!   Chris Gaskill, 7/89  Added early Tobiska model
!   Stan Solomon,  8/89  Corrections to above
!   Stan Solomon,  1/90  Tobiska SERF2; added W & R spectra
!   Stan Solomon,  6/91  Tobiska EUV 91; Hntggr Ly-B, Fe XVI scaling
!   Stan Solomon,  2/92  Updated Tobiska EUV91; corrected SME FUV
!   Scott Bailey, 12/93  Initial one-nm bins version
!   Stan Solomon,  6/04  Added EUVAC option, cleaned up artifacts
!   Stan Solomon,  9/04  Added ability to specify input data file
!   Stan Solomon,  3/05  Changed all to photon units
!   Stan Solomon,  1/15  Updated for f90; only read file on first call
!
! Calling parameters:
! ISCALE   =0 for Hinteregger contrast ratio method
!          =1 for EUVAC
!          =2 for user-supplied data
! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
! F107A    81-day centered average 10.7 cm flux
! HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
! FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
! HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
! HEIEW    He I 10830A equivalent width (mAngstroms) (obsolete)
! XUVFAC   factor for scaling flux 18-250A or 18-50A (optional)

! Returned parameters:
! WAVE1    longwave bound of spectral intervals (Angstroms)
! WAVE2    shortwave bound of intervals
! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
!
! Other definitions:
! LMAX     dimension of flux and scaling arrays, currently = 123
! WAVEL    = WAVE1
! WAVES    = WAVE2
! RFLUX    low solar activity flux
! XFLUX    high solar activity flux
! SCALE1   scaling factors for H LyB-keyed chromospheric emissions
! SCALE2   scaling factors for FeXVI-keyed coronal emissions
! B1       fit coefficients for H LyB
! B2       fit coefficients for FeXVI
! R1       enhancement ratio for H LyB
! R2       enhancement ratio for FeXVI
! P107     average of F107 and F107A
! A        scaling factor for EUVAC model
!
!
      SUBROUTINE SSFLUX (ISCALE, F107, F107A, HLYBR, FEXVIR, HLYA,
     >                   HEIEW, XUVFAC, WAVE1, WAVE2, SFLUX)
!      use cglow,only: lmax
      implicit none
      include 'cglow.h'

! Args:
      integer,intent(in) :: iscale
      real,intent(in)    :: f107, f107a,HLYBR, FEXVIR,HLYA,XUVFAC,heiew
      real,intent(out),dimension(Lmax)   :: wave1,wave2,sflux
! Local:
      real WAVEL(LMAX), WAVES(LMAX), RFLUX(LMAX), UFLUX(LMAX),
     >          SCALE1(LMAX), SCALE2(LMAX), A(LMAX),p107,r1,r2
      integer l
      real,parameter :: epsil=1.0E-6
!      integer :: islast=-1  !this didn't work right in f2py


! regression coefficients which reduce to solar min. spectrum:
      real :: B1(NMAJ),B2(NMAJ)
      DATA B1/1.0, 0.0138, 0.005/
      DATA B2/1.0, 0.59425, 0.3811/

! 'best fit' regression coefficients, commented out, for reference:
!     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/


! Hinteregger contrast ratio method:

      IF (iscale .eq. 0) then
!        if (islast .ne. iscale) then
          open(unit=1,file='ssflux_hint.dat',status='old')
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),rflux(l),scale1(l),scale2(l)
          enddo
          close(unit=1)
!         endif

        IF (HLYBR .GT. EPSIL) THEN
          R1 = HLYBR
        ELSE
          R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
        ENDIF
        IF (FEXVIR .GT. EPSIL) THEN
          R2 = FEXVIR
        ELSE
          R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
        ENDIF

        do l=1,lmax
          SFLUX(L) = RFLUX(L) + (R1-1.)*SCALE1(L) + (R2-1.)*SCALE2(L)
          IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
          IF (XUVFAC .GT. EPSIL .AND.
     >        WAVEL(L).LT.251.0 .AND. WAVES(L).GT.17.0)
     >        SFLUX(L)=SFLUX(L)*XUVFAC
        enddo
      endif

! EUVAC Method:

      IF (iscale .eq. 1) then
!        if (islast .ne. iscale) then
          open(unit=1,file='ssflux_euvac.dat',status='old')
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),rflux(l),a(l)
          enddo
          close(unit=1)
!        endif

      P107 = (F107+F107A)/2.
        do l=1,lmax
          SFLUX(L) = RFLUX(L) * (1. + A(L)*(P107-80.))
          IF (SFLUX(L) .LT. 0.8*RFLUX(L)) SFLUX(L) = 0.8*RFLUX(L)
          IF (XUVFAC .GT. EPSIL .AND.
     >        WAVEL(L).LT.51.0 .AND. WAVES(L).GT.17.0)
     >        SFLUX(L)=SFLUX(L)*XUVFAC
        enddo
      ENDIF

! User-supplied data:

      if (iscale .eq. 2) then
!        if (islast .ne. iscale) then
          open(unit=1,file='ssflux_user.dat',status='old')
          read(1,*)
          do l=lmax,1,-1
            read(1,*) waves(l),wavel(l),uflux(l)
          enddo
          close(unit=1)
!        endif
        do l=1,lmax
          sflux(l)=uflux(l)
        enddo
      endif

! Fill wavelength arrays, substitute in H Lyman-alpha if provided:

      do l=1,lmax
        WAVE1(L) = WAVEL(L)
        WAVE2(L) = WAVES(L)
        IF (HLYA .GT. EPSIL .AND.
     >      WAVEL(L).LT.1221. .AND. WAVES(L).GT.1209.)
     >      SFLUX(L) = HLYA
      enddo

!      islast=iscale

      End SUBROUTINE SSFLUX
