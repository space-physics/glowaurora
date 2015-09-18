! Subroutine MAXT
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Stan Solomon, 11/89, 9/91, 1/94, 3/05
!
! Generates Maxwellian electron spectra with, optionally, a low energy tail
! of the form used by Meier et al., JGR 94, 13541, 1989.
!
! Supplied by calling routine:
!     EFLUX  total energy flux in erg cm-2 s-1
!     EZER   characteristic energy in eV
!     ENER   energy grid in eV
!     dE   energy bin width in eV
!     NBINS  number of energy bins (dimension of ENER, dE, and PHI)
!     ITAIL  1 = Maxwellian with low-energy tail, 0 = regular Maxwellian
!     FMONO  additional monoenergetic energy flux in erg cm-2 s-1
!     EMONO  characteristic enerngy of FMONO in eV
!
! Returned by subroutine:
!     MAXT   Hemispherical flux in cm-2 s-1 eV-1
!
      Subroutine maxt(EFLUX, EZER, ENER, dE, ITAIL, FMONO, EMONO, phi)
!      use cglow, only: nbins
      implicit none
      include 'cglow.h'
!Args:
      Real,Intent(Out) :: phi(NBINS)
      Real, Intent(In)  :: EFLUX, EZER, ENER(NBINS), dE(NBINS),
     & FMONO,EMONO
      Integer,Intent(In):: ITAIL
!Local:
      Real B,TE, PHIMAX, ERAT
      Integer K

      TE = 0.
C
      IF (EZER < 500.) THEN
        B = 0.8*EZER
      ELSE
        B = 0.1*EZER + 350.
      ENDIF
C
      PHIMAX = EXP(-1.)
C
      DO 300 K=1,NBINS
        ERAT = ENER(K) / EZER
        IF (ERAT > 60.) ERAT = 60.
        PHI(K) = ERAT * EXP(-ERAT)
        IF (ITAIL > 0)
     >    PHI(K) = PHI(K) + 0.4*PHIMAX*(EZER/ENER(K))*EXP(-ENER(K)/B)
        TE = TE + PHI(K) * dE(K) * ENER(K) * 1.6022E-12
  300 CONTINUE
C
      DO 400 K=1,NBINS
      PHI(K) = PHI(K) * EFLUX / TE
  400 CONTINUE
c
c
c
      if (fmono > 0.) then
      do 500 k=1,nbins
      if (emono > ener(k)-dE(k)/2..and.emono < ener(k)+dE(k)/2.)
     >  phi(k)=phi(k)+fmono/(1.6022E-12*dE(k)*ener(k))
  500 continue
      endif
      END Subroutine maxt
