! Subroutine SNOEMINT gets NO estimate from the NOEM emperical model and
! INTerpolates it onto an altitude grid.  Extrapolation is done above 150
! km assuming a scale height approximation, and below 100 km
! assuming a constant profile.
!
! Stan Solomon, 12/14
!
! Input:
!   IDATE  Date in yyddd or yyyyddd format
!   GLAT   Geographic latitude in degrees
!   GLONG  Geographic longitude in degrees
!   F107   10.7 cm radio flux index
!   AP     Ap index
!   JMAX   Number of points in altitude grid
!   Z      Altitude grid in km
!   ZTN    Temperature at Z in K
! Output:
!   ZNO    Nitric oxide density at Z in cm-3
!
!
module SNOEMINTmod
    use ccglow
    use snoemmod
    implicit none
    private
    public :: snoemint

contains
      SUBROUTINE SNOEMINT(IDATE,GLAT,GLONG,F107,AP,Z,ZTN,ZNO)

      integer, intent(in) :: IDATE
      real(sp),intent(in) :: GLAT,GLONG,F107,AP,Z(JMAX),ZTN(JMAX)
      real(sp),intent(out):: ZNO(JMAX)

      Real(sp) :: ZG(16),XMLATNO(33), ZMNO(33,16), ZMNOI(16),rat,xkp,xmlat,xmlong,H
      integer :: iday,j,klat1,klat2,kz1,kz2

! Find magnetic latitude:
      CALL GEOMAG(0,GLONG,GLAT,XMLONG,XMLAT)

! Get zonal mean NO profiles:
      IDAY=IDATE-IDATE/1000*1000
      XKP=1.75*LOG(0.4*AP)
      CALL SNOEM(IDAY,XKP,F107,ZG,XMLATNO,ZMNO)

! Interpolate altitude profile at magnetic latitude:
      KLAT1=INT(XMLAT+80.)/5+1
      KLAT2=KLAT1+1
      IF (KLAT1 .LT. 1) KLAT1=1
      IF (KLAT1 .GT. 33) KLAT1=33
      IF (KLAT2 .LT. 1) KLAT1=1
      IF (KLAT2 .GT. 33) KLAT2=33
      RAT=XMLAT/5.-INT(XMLAT)/5
!
      DO J=1,16
        ZMNOI(J) = LOG(ZMNO(KLAT1,J)*(1.-RAT)+ZMNO(KLAT2,J)*RAT)
      END DO
!
      H=0.03*ZTN(JMAX)
      DO J=1,JMAX
        IF (Z(J) .LE. 100.) ZNO(J)=EXP(ZMNOI(16))
        IF (Z(J) .GT. 100. .AND. Z(J) .LE. 150.) THEN
          KZ2=INT((150.-Z(J))*.3)+1
          KZ1=KZ2+1
          ZNO(J)=EXP ( ZMNOI(KZ1) + (ZMNOI(KZ2)-ZMNOI(KZ1)) &
                                * (Z(J)-ZG(KZ1)) / (ZG(KZ2)-ZG(KZ1)) )
        ENDIF
        IF (Z(J) .GT. 150.) ZNO(J)=EXP(ZMNOI(1)+(150.-Z(J))/H)
      END DO
!
      END subroutine snoemint
end module SNOEMINTmod
