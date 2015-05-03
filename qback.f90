! Subroutine QBACK
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Stan Solomon, 11/88, 11/92
! Comment updated 3/05
!
! Estimates background ("nighttime") ionization rates.
! Intended to conform to the TIGCM background (Roble et al., 1989).
! Three fluxes, from the stellar background and multiple scattering
! of solar atomic hydrogen emissions in the geocorona, are considered:
! The sum of ionizing flux in the 100-1000 A range (FIONT)
! H Lyman beta (FLYBT)
! H Ly alpha (FLYAT)
!
module qbackmod
    use ccglow
    implicit none
    private
    public :: qback

  real(sp),parameter :: FIONT=5.0E7, FLYBT=1.0E7, FLYAT=1.0E9, SIGIO=1.0E-17, &
  SIGIO2=2.0E-17, SIGIN2=2.0E-17, SLBAO2=1.6E-18, SLBIO2=1.0E-18, &
  SLAAO2=1.0E-20, SLAINO=2.0E-18
contains
  SUBROUTINE QBACK (ZMAJ, ZNO, ZVCD, PHOTOI, PHONO)

  real(sp), intent(in) :: ZMAJ(NMAJ,Jmax),ZNO(:),ZVCD(NMAJ,Jmax)
  !PHOTI and PHONO must be INOUT
  real(sp),intent(INout)::PHOTOI(NST,NMAJ,Jmax), PHONO(NST,Jmax)

  real(sp), dimension(jmax) :: taui, taulya,taulyb, fion
!
! Calculate ionization rates at each altitude:
!
    TAUI = 2.*(SIGIO*ZVCD(1,:)+SIGIO2*ZVCD(2,:)+SIGIN2*ZVCD(3,:))
    Where (TAUI > 60.) TAUI = 60.
    TAULYB = 2.*SLBAO2*ZVCD(2,:)
    Where (TAULYB > 60.) TAULYB = 60.
    TAULYA = 2.*SLAAO2*ZVCD(2,:)
    Where (TAULYA > 60.) TAULYB = 60.
    FION = FIONT * EXP(-TAUI)
    PHOTOI(1,1,:) = PHOTOI(1,1,:) + FION * ZMAJ(1,:) * SIGIO
    PHOTOI(1,2,:) = PHOTOI(1,2,:) + FION * ZMAJ(2,:) * SIGIO2 &
                   + FLYBT * EXP(-TAULYB) * ZMAJ(2,:) * SLBIO2
    PHOTOI(1,3,:) = PHOTOI(1,3,:) + FION * ZMAJ(3,:) * SIGIN2
    PHONO(1,:) = PHONO(1,:) + FLYAT * EXP(-TAULYA) * ZNO(:) * SLAINO
!
      END SUBROUTINE QBACK
end module qbackmod
