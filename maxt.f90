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
module MAXT
 use machprec
 public :: phi0
contains
 Subroutine phi0(EFLUX, EZER, ENER, dE, NBINS,ITAIL, FMONO, EMONO, phi)
    implicit none
    Real(sp),Intent(Out) :: phi(NBINS)

    Real(sp), Intent(In)  :: EFLUX, EZER, ENER(NBINS), dE(NBINS), FMONO,EMONO
    Integer,Intent(In):: NBINS,ITAIL

    Real(sp) :: B,TE, PHIMAX, ERAT
    Integer :: K

    !
    TE = 0.
    !
    IF (EZER < 500.) THEN
    B = 0.8*EZER
    ELSE
    B = 0.1*EZER + 350.
    ENDIF
    !
    PHIMAX = EXP(-1.)
    !
    DO K=1,NBINS
        ERAT = ENER(K) / EZER
        IF (ERAT > 60.) ERAT = 60.
        phi(K) = ERAT * EXP(-ERAT)
        IF (ITAIL > 0) then
          phi(K) = phi(K) + 0.4*PHIMAX*(EZER/ENER(K))*EXP(-ENER(K)/B)
        End If
        TE = TE + phi(K) * dE(K) * ENER(K) * 1.6022E-12
    End Do
    !
    phi = phi * EFLUX / TE
    !
    !
    if (fmono > 0.) then
        Do k=1,nbins
          if (emono > ener(k)-dE(k)/2. .and. emono < ener(k)+dE(k)/2.) Then
             phi(k)=phi(k)+fmono/(1.6022E-12*dE(k)*ener(k))
          End If
        End Do
    End If
    !

 END Subroutine phi0
end module MAXT
