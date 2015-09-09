! Subroutine SOLZEN
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Stan Solomon, 1988
! Temporary Y2K fix-up, SCS, 2005.
!
! Returns Solar Zenith Angle SZA in degrees for specified date in form
! yyddd, universal time in seconds, geographic latitude and longitude
! in degrees.

      SUBROUTINE SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
!      use cglow,only: pi
      implicit none
      include 'cglow.h'
! Args:
      integer,intent(in) :: idate
      real,intent(in) :: UT,glat,glong
      real,intent(out):: sza
! Local:  
      real sdec,srasn,gst,rlat,rh,cossza,rlong

      RLAT = GLAT * PI/180.
      RLONG = GLONG * PI/180.
      CALL SUNCOR (IDATE, UT, SDEC, SRASN, GST)
      RH = SRASN - (GST+RLONG)
      COSSZA = SIN(SDEC)*SIN(RLAT) + COS(SDEC)*COS(RLAT)*COS(RH)
      SZA = ACOS(COSSZA) * 180./PI
      END SUBROUTINE SOLZEN
!
! Subroutine SUNCOR returns the declination SDEC and right ascension
! SRASN of the sun in GEI coordinates, radians, for a given date IDATE
! in yyddd format and universal time UT in seconds.  Greenwich Sidereal
! Time GST in radians is also returned.  Reference:  C.T. Russell,
! Geophysical Coordinate Transforms.
!
      SUBROUTINE SUNCOR (IDATE, UT, SDEC, SRASN, GST)
!      use cglow,only: pi
      implicit none
      include 'cglow.h'
! Args:
      integer,intent(in) :: idate
      real,intent(in) :: UT
      real,intent(out):: SDEC,SRASN,GST
! Local:
      real fday,vl,slp,slong,sind,obliq,GG,DJ,cosd,T
      integer iyr,iday

      FDAY=UT/86400.
      IYR=IDATE/1000
      IDAY=IDATE-IYR*1000
!
! Temporary Y2K fix-up:
! Should work with either yyddd or yyyyddd format from 1950 to 2050.
! Note deteriorating accuracy after ~2050 anyway.
! Won't work after 2100 due to lack of a leap year.
      IF (IYR .GE. 1900) IYR=IYR-1900
      IF (IYR .LT. 50) IYR=IYR+100

      DJ=365*IYR+(IYR-1)/4+IDAY+FDAY-0.5
      T=DJ/36525.
      VL= MOD(279.696678+.9856473354*DJ,360.)
      GST=MOD(279.696678+.9856473354*DJ+360.*FDAY+180.,360.) * PI/180.
      GG=  MOD(358.475845+.985600267*DJ,360.) * PI/180.
      SLONG=VL+(1.91946-.004789*T)*SIN(GG)+.020094*SIN(2.*GG)
      OBLIQ=(23.45229-0.0130125*T) *PI/180.
      SLP=(SLONG-.005686) * PI/180.
      SIND=SIN(OBLIQ)*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SDEC=ATAN(SIND/COSD)
      SRASN=pi-ATAN2(1./TAN(OBLIQ)*SIND/COSD,-COS(SLP)/COSD)
      END Subroutine Suncor
