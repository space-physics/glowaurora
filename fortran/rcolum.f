! Subroutine RCOLUM
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Stan Solomon, 1988, 1991
!
! Calculates the column density ZCOL for each species ZMAJ above height
! ZZ at zenith angle CHI.  Uses a Chapman function fit [Smith and Smith,
! JGR 77, 3592, 1972].  If CHI is less than 90 degrees, column
! densities are calculated directly; if CHI is greater than 90 degrees
! the column density at grazing height for 90 degrees is calculated and
! doubled and the column density above ZZ(J) is subtracted; if CHI is
! such that the grazing height is less than the radius of the earth the
! column densities are set to 'infinity', i.e., 1.0E30.  Densities
! supplied in array ZMAJ are used in the calculation except where
! grazing height is below the lowest level specified, in this case
! values are interpolated logarithmically from a US Standard Atmosphere
! at sea level, the tropopause, the stratopause, and the mesopause.
!
      SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, ZCOL, ZVCD)
!      use cglow,only: nmaj,jmax,pi,Re
      implicit none
      include 'cglow.h'
!Args:
      real, intent(in) :: chi,zz(jmax),ZMAJ(NMAJ,JMAX),TN(JMAX)
      real, intent(out) :: ZVCD(NMAJ,JMAX),ZCOL(NMAJ,JMAX)
!Local:
      integer,PARAMETER :: NM=3,NU=4

      real  ZCG(NM),ghrg,ghz,tng
      real  ZUS(nu), TNUS(nu),zcus(nm,nu)
      integer i,j, jg

      real :: chap !function

      DATA ZUS/0.0, 1.5E6, 5.E6, 9.E6/
      DATA TNUS/288., 217., 271., 187./
      DATA ZCUS/8.00E17, 4.54E24, 1.69E25,
     >          8.00E17, 5.46E23, 2.03E24,
     >          8.00E17, 3.63E21, 1.35E22,
     >          7.80E17, 8.48E18, 3.16E19/

      CALL VCD (ZZ, ZMAJ, ZVCD)

      IF (CHI .GE. 2.) THEN 
        DO 40 I=1,NMAJ
        DO 40 J=1,JMAX
        ZCOL(I,J) = 1.0E30
   40   CONTINUE
        RETURN
      ENDIF
C
      IF (CHI .LE. PI/2.) THEN
        DO 60 I=1,NMAJ
        DO 60 J=1,JMAX
        ZCOL(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
   60   CONTINUE
      ELSE
        DO 220 J=1,JMAX
        GHRG=(RE+ZZ(J))*SIN(CHI) 
        GHZ=GHRG-RE 
        IF (GHZ .LE. 0.) THEN
          DO 80 I=1,NMAJ
          ZCOL(I,J) = 1.0E30
   80     CONTINUE
          GOTO 220
        ENDIF
        IF (GHZ .GE. ZZ(1)) THEN
          DO 100 JG=1,J-1
            IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
  100     CONTINUE
  120     TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
          DO 140 I=1,NMAJ
          ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) **
     >                      ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
  140     CONTINUE
        ELSE
          DO 160 JG=1,3
          IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
  160     CONTINUE
  180     TNG = TNUS(JG)
     >        + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
          DO 200 I=1,NMAJ
          ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) **
     >                      ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
  200     CONTINUE
        ENDIF
        DO 210 I=1,NMAJ
        ZCOL(I,J) = 2. * ZCG(I) * CHAP(PI/2.,GHZ,TNG,I)
     >              - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
  210   CONTINUE
  220   CONTINUE
      ENDIF
      END SUBROUTINE RCOLUM


      real FUNCTION CHAP (CHI, Z, T, I)
!      use cglow,only: nmaj,pi,Re,G
      implicit none
      include 'cglow.h'
!Args:
      real,intent(in) :: chi,z,t
      integer,intent(in) :: I
!Local:
      real AM(nmaj)
      real gr,hn,hg,hf,sqhf,SPERFC 

      DATA AM/16., 32., 28./

      GR=G*(RE/(RE+Z))**2 
      HN=1.38E-16*T/(AM(I)*1.662E-24*GR)
      HG=(RE+Z)/HN 
      HF=0.5*HG*(COS(CHI)**2) 
      SQHF=SQRT(HF) 
      CHAP=SQRT(0.5*PI*HG)*SPERFC(SQHF) 
      End Function CHAP


      pure real FUNCTION SPERFC(DUMMY) 
      implicit none
      real,intent(in) :: dummy
      IF (DUMMY .LE. 8.) THEN
        SPERFC = (1.0606963+0.55643831*DUMMY) /
     >           (1.0619896+1.7245609*DUMMY+DUMMY*DUMMY)
      ELSE
        SPERFC=0.56498823/(0.06651874+DUMMY) 
      ENDIF 

      END function sperfc


      SUBROUTINE VCD(ZZ,ZMAJ,ZVCD)
!      use cglow,only: nmaj,jmax
      implicit none
      include 'cglow.h'
! Args:
      real,intent(in) :: zz(jmax),zmaj(NMAJ,JMAX)
      real,intent(out) :: ZVCD(NMAJ,JMAX)
! Local:
      real rat
      integer i, j

      DO 200 I=1,NMAJ
      ZVCD(I,JMAX) =   ZMAJ(I,JMAX)
     >               * (ZZ(JMAX)-ZZ(JMAX-1))
     >               / LOG(ZMAJ(I,JMAX-1)/ZMAJ(I,JMAX))
      DO 200 J=JMAX-1,1,-1
      RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
      ZVCD(I,J) =   ZVCD(I,J+1)
     >            + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / LOG(RAT) * (1.-RAT)
  200 CONTINUE
      END Subroutine VCD
