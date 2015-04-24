C Subroutine EGRID sets up electron energy grid
C
C Stan Solomon, 1/92
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
      SUBROUTINE EGRID (ENER, DEL, NBINS)
      Implicit None
      Integer N
      Integer, Intent(In)  :: NBINS
      Real, Intent(Out) :: ENER(NBINS), DEL(NBINS)

      DO N=1,NBINS
        IF (N .LE. 21) THEN
          ENER(N) = 0.5 * FLOAT(N)
        ELSE
          ENER(N) = EXP (0.05 * FLOAT(N+26))
        ENDIF
      End Do

      DEL(1) = 0.5

      DO N=2,NBINS
        DEL(N) = ENER(N)-ENER(N-1)
      End Do

      DO N=1,NBINS
        ENER(N) = ENER(N) - DEL(N)/2.0
      End Do

      END
