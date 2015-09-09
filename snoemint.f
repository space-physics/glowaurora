C Subroutine SNOEMINT gets NO estimate from the NOEM emperical model and
C INTerpolates it onto an altitude grid.  Extrapolation is done above 150
C km assuming a scale height approximation, and below 100 km
C assuming a constant profile.
C
C Stan Solomon, 12/14
C
C Input:
C   IDATE  Date in yyddd or yyyyddd format
C   GLAT   Geographic latitude in degrees
C   GLONG  Geographic longitude in degrees
C   F107   10.7 cm radio flux index
C   AP     Ap index
C   JMAX   Number of points in altitude grid
C   Z      Altitude grid in km
C   ZTN    Temperature at Z in K
C Output:
C   ZNO    Nitric oxide density at Z in cm-3
C
C
      SUBROUTINE SNOEMINT(IDATE,GLAT,GLONG,F107,AP,Z,ZTN,ZNO)
!      use cglow,only: jmax
      implicit none
      include 'cglow.h'
      
      integer, intent(in) :: IDATE
      real,intent(in) :: GLAT,GLONG,F107,AP,Z(JMAX),ZTN(JMAX)
      real,intent(out):: ZNO(JMAX)
      
      Real ZG(16),XMLATNO(33), ZMNO(33,16), ZMNOI(16),rat,xkp,
     & xmlat,xmlong,h
      integer iday,j,klat1,klat2,kz1,kz2
      

C Find magnetic latitude:

      CALL GEOMAG(0,GLONG,GLAT,XMLONG,XMLAT)
C
C
C Get zonal mean NO profiles:
C
      IDAY=IDATE-IDATE/1000*1000
      XKP=1.75*LOG(0.4*AP)
      CALL SNOEM(IDAY,XKP,F107,ZG,XMLATNO,ZMNO)
C
C
C Interpolate altitude profile at magnetic latitude:
C
      KLAT1=INT(XMLAT+80.)/5+1
      KLAT2=KLAT1+1
      IF (KLAT1 .LT. 1) KLAT1=1
      IF (KLAT1 .GT. 33) KLAT1=33
      IF (KLAT2 .LT. 1) KLAT1=1
      IF (KLAT2 .GT. 33) KLAT2=33
      RAT=XMLAT/5.-INT(XMLAT)/5
C
      DO J=1,16
        ZMNOI(J) = LOG(ZMNO(KLAT1,J)*(1.-RAT)+ZMNO(KLAT2,J)*RAT)
      END DO

      H=0.03*ZTN(JMAX)
      DO J=1,JMAX
        IF (Z(J) .LE. 100.) ZNO(J)=EXP(ZMNOI(16))
        IF (Z(J) .GT. 100. .AND. Z(J) .LE. 150.) THEN
          KZ2=INT((150.-Z(J))*.3)+1
          KZ1=KZ2+1
          ZNO(J)=EXP ( ZMNOI(KZ1) + (ZMNOI(KZ2)-ZMNOI(KZ1))
     >                            * (Z(J)-ZG(KZ1)) / (ZG(KZ2)-ZG(KZ1)) )
        ENDIF
        IF (Z(J) .GT. 150.) ZNO(J)=EXP(ZMNOI(1)+(150.-Z(J))/H)
      END DO
      END SUBROUTINE SNOEMINT
