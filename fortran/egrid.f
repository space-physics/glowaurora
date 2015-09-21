! Subroutine EGRID sets up electron energy grid
!
! Stan Solomon, 1/92
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
      SUBROUTINE EGRID (ENER, DEL)
!      use cglow,only: nbins
      implicit none
      include 'cglow.h'
! Args:
      Real, Intent(Out) :: ENER(Nbins), DEL(Nbins)
! Local:
      Integer N

      DO N=1,NBINS
        IF (N .LE. 21) THEN
          ENER(N) = 0.5 * REAL(N)
        ELSE
          ENER(N) = EXP (0.05 * REAL(N+26))
        ENDIF
      end do

      DEL(1) = 0.5
      DO N=2,NBINS
        DEL(N) = ENER(N)-ENER(N-1)
      End Do

      DO N=1,NBINS
        ENER(N) = ENER(N) - DEL(N)/2.0
      End Do

      End Subroutine Egrid
