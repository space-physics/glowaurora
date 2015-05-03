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
module rcolummod
  use machprec
  implicit none

  private
   Real(sp),parameter :: RE=6.37E8, G=978.1, kerg=1.38E-16,amug=1.662E-24
   integer, parameter:: NM=3,NU=4
   Real(sp),parameter:: ZCUS(NM,NU)=reshape( &
                                   [8.00E17, 4.54E24, 1.69E25, 8.00E17, &
                                   5.46E23, 2.03E24, 8.00E17, 3.63E21, &
                                   1.35E22, 7.80E17, 8.48E18, 3.16E19],&
                                   [NM,NU])
   Real(sp), parameter,dimension(nu) :: ZUS=[0., 1.5E6, 5.E6, 9.E6], &
                                        TNUS=[288., 217., 271., 187.]
   Real(sp), parameter :: AM(nm)=[16., 32., 28.]

   public :: rcolum, chap,vcd !f2py needs all subr/func public

contains
  SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, ZCOL, ZVCD, NMAJ)
  include 'glow.h'
  integer,intent(in) :: nmaj
  real(sp), intent(in) :: chi,zz(jmax),TN(JMAX)
  real,intent(out) :: ZVCD(NMAJ,JMAX),ZCOL(NMAJ,JMAX)


!
  real(sp) :: ZMAJ(NMAJ,JMAX), ZCG(NM), ghrg,ghz,tng
  integer i,j, jg
!


!
!
  CALL VCD (ZZ, ZMAJ, ZVCD, JMAX, NMAJ)
!
  IF (CHI .GE. 2.) THEN
    ZCOL = 1.0E30
    RETURN
  ENDIF
!
  IF (CHI .LE. PI/2.) THEN
    DO 60 J=1,JMAX
    ZCOL(:,J) = ZVCD(:,J) * CHAP(CHI,ZZ(J),TN(J))
60   CONTINUE
  ELSE
    DO 220 J=1,JMAX
    GHRG=(RE+ZZ(J))*SIN(CHI)
    GHZ=GHRG-RE
    IF (GHZ .LE. 0.) THEN
      DO I=1,NMAJ
          ZCOL(I,J) = 1.0E30
      End Do
      GOTO 220
    ENDIF
    IF (GHZ .GE. ZZ(1)) THEN
      DO JG=1,J-1
        IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) EXIT
      End DO
      TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
      DO I=1,NMAJ
      ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG))**((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
      End Do
    ELSE
      DO JG=1,3
        IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) Exit
      End Do
      TNG = TNUS(JG) + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))

      ZCG = ZCUS(:,JG) * (ZCUS(:,JG+1) / ZCUS(:,JG))**((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))

    ENDIF

    ZCOL(:,J) = 2. * ZCG * CHAP(PI/2.,GHZ,TNG) - ZVCD(:,J) * CHAP(CHI,ZZ(J),TN(J))

220   CONTINUE
  ENDIF
!
  END SUBROUTINE RCOLUM
!
!
  FUNCTION CHAP (CHI, Z, T)
      real(sp),intent(in) :: chi,z,t
      real(sp),allocatable,dimension(:) :: gr,hn,hg,hf,sqhf,SPERFC,chap

      GR=G*(RE/(RE+Z))**2
      HN=kerg*T/(AM *amug*GR)
      HG=(RE+Z)/HN
      HF=0.5*HG*(COS(CHI)**2)
      SQHF=SQRT(HF)
      Where(SQHF <= 8.)
       SPERFC = (1.0606963+0.55643831*SQHF) / (1.0619896+1.7245609*SQHF+SQHF**2)
      ELSEWhere
       SPERFC=0.56498823/(0.06651874+SQHF)
      END Where
      CHAP=SQRT(0.5*PI*HG)*SPERFC
  END Function Chap
!
!
!
  SUBROUTINE VCD(ZZ,ZMAJ,ZVCD,JMAX,NMAJ)
  integer,intent(in) :: nmaj,JMAX
  real(sp),intent(in) :: zz(:),zmaj(NMAJ,JMAX)
  real(sp),intent(out) :: ZVCD(NMAJ,JMAX)
  real(sp) rat
  integer i, j

  DO I=1,NMAJ
  ZVCD(I,JMAX) =   ZMAJ(I,JMAX) * (ZZ(JMAX)-ZZ(JMAX-1)) / LOG(ZMAJ(I,JMAX-1)/ZMAJ(I,JMAX))
      DO J=JMAX-1,1,-1
          RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
          ZVCD(I,J) =   ZVCD(I,J+1) + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / LOG(RAT) * (1.-RAT)
      End Do
  End Do
  END Subroutine VCD
end module rcolummod
