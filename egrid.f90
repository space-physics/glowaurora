! Subroutine EGRID sets up electron energy grid
!
! Stan Solomon, 1/92
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
module EnergyGrid
  use machprec
  public :: EGRID
contains
  SUBROUTINE EGRID (ENER, DEL,NBINS)
      Implicit None
      Integer N
      Integer,Intent(In) :: Nbins
      Real(sp), Intent(Out) :: ENER(Nbins), DEL(Nbins)


      DO N=1,NBINS
        IF (N .LE. 21) THEN
          ENER(N) = 0.5 * FLOAT(N)
        ELSE
          ENER(N) = EXP (0.05 * FLOAT(N+26))
        ENDIF
      End Do

      DEL(1) = 0.5
      DEL(2:nbins) = ENER(2:nbins)-ENER(1:nbins-1)

      ENER = ENER - DEL/2.0

  END Subroutine EGRID
end module EnergyGrid
