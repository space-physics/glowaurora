! Subroutine EGRID sets up electron energy grid
!
! Stan Solomon, 1/92
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
module energygrid
  use ccglow
  Implicit None
  private
  public :: EGRID
contains
  SUBROUTINE EGRID (ENER, DEL,NBINS)

      Integer,Intent(In) :: Nbins
      Real(sp), Intent(Out) :: ENER(Nbins), DEL(Nbins)
      Real :: tmp(Nbins)
      Integer :: i,N(NBINS)
       N = [(i,i=1,nbins)] !must be on a separate line

tmp = 0.05 * real(N+26)
      !print*,maxexponent(ener)
        Where (N <= 21)
          ENER = 0.5 * REAL(N)
        ELSEwhere

          ENER = EXP (tmp)
        END where

      DEL(1) = 0.5
      DEL(2:nbins) = ENER(2:nbins)-ENER(1:nbins-1)

      ENER = ENER - DEL/2.0
    END Subroutine EGRID
end module EnergyGrid
