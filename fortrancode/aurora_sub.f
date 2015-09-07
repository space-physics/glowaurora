C Example driver program for GLOW subroutine package - aurora version
C
C Stan Solomon, 4/05, 12/14
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C
C For more information see the file glow.txt
C
C For definitions of common block /CGLOW/, see subroutine GLOW
C
C Other definitions:
C F107P   Solar 10.7 cm flux for previous day
C AP      Ap index of geomagnetic activity
C Z       altitude array, km
C XNO     default values for NO concentration, cm-3
C
C Array dimensions:
C JMAX    number of altitude levels
C NBINS   number of energetic electron energy bins
C LMAX    number of wavelength intervals for solar flux
C NMAJ    number of major species
C NEX     number of ionized/excited species
C NW      number of airglow emission wavelengths
C NC      number of component production terms for each emission
C NST     number of states produced by photoionization/dissociation
C NEI     number of states produced by electron impact
C NF      number of types of auroral fluxes
C
      SUBROUTINE AURORA(Z,Pyion,Pyecalc,Pypi,Pysi,Pyisr,
     &                  Pyidate, Pyut, Pyglat, Pyglong, Pyf107a, Pyf107,
     &                  Pyf107p, Pyap,PyPhitop)

      use cglow,only: nbins,jmax,nmaj

      Integer, Intent(In) :: Pyidate
      Real,Intent(In) :: Pyut, Pyglat, Pyglong, Pyf107a, Pyf107,
     &                  Pyf107p, Pyap, PyPhitop(NBINS,nmaj)
      Real, Dimension(JMAX),Intent(Out)    :: Pyecalc,Pypi,Pysi
      Real, Intent(Out)  :: Pyion(JMAX,11), Pyisr(JMAX,nmaj)

C
      COMMON /CGLOW/
     >    IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM,
     >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
     >    ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
     >    ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
     >    ZTN(JMAX), ZTI(JMAX), ZTE(JMAX),
     >    PHITOP(NBINS), EFLUX(NF), EZERO(NF),
     >    SZA, DIP, EFRAC, IERR,
     >    ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX),
     >    WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
     >    ENER(NBINS), DEL(NBINS),
     >    PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
     >    PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX),
     >    QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX),
     >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
     >    EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX),
     >    ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW)
C
      COMMON /CXSECT/ SIGS(NMAJ,NBINS), PE(NMAJ,NBINS), PIN(NMAJ,NBINS),
     >                SIGA(NMAJ,NBINS,NBINS), SEC(NMAJ,NBINS,NBINS),
     >                SIGEX(NEI,NMAJ,NBINS), SIGIX(NEI,NMAJ,NBINS),
     >                IIMAXX(NBINS)
C
      DIMENSION Z(JMAX), D(8), T(2), SW(25),
     >          OUTF(11,JMAX), OARR(30), TPI(NMAJ)
C
      LOGICAL JF(12)
C
      DATA SW/25*1./
C     Hack alert, didn't get fancier due to newer GLOW version coming soon enough
      idate=Pyidate; ut=Pyut; glat=Pyglat; glong=Pyglong;
      f107a=Pyf107a; f107=Pyf107; f107p=Pyf107p; ap=Pyap
      ENER = PyPhitop(:,1); DEL=PyPhitop(:,2); PHITOP=PyPhitop(:,3)
C
C Obtain input parameters:
C
C     read (5,*) idate, ut, glat, glong, f107a, f107, f107p, ap, ef, ec
C
C
C Set other parameters and switches:
C
      JLOCAL = 0
      KCHEM=4
      ISCALE = 1
      XUVFAC = 3.
      HLYBR = 0.
      FEXVIR = 0.
      HLYA = 0.
      HEIEW = 0.
      ITAIL = 0
      FMONO = 0.
      EMONO = 0.
C
C Calculate local solar time:
C
      STL = (UT/240.+GLONG) / 15.
      IF (STL .LT. 0.) STL = STL + 24.
      IF (STL .GT. 24.) STL = STL - 24.
C
C
C Call MSIS-2K to get neutral densities and temperature:
C
        CALL TSELEC(SW)

        DO J=1,JMAX
          CALL GTD7(IDATE,UT,Z(J),GLAT,GLONG,STL,F107A,F107P,AP,48,D,T)
          ZO(J) = D(2)
         ZN2(J) = D(3)
         ZO2(J) = D(4)
          ZRHO(J) = D(6)
         ZNS(J) = D(8)
          ZTN(J) = T(2)
       END DO
C
C
C Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical
C Model (NOEM)
C
      CALL SNOEMINT(IDATE,GLAT,GLONG,F107,AP,Z,ZTN,ZNO)
C
C
C Call International Reference Ionosphere-1990 subroutine to get
C electron density and temperature and ion temperature:
C
C NOTE: the directory specified in the call to IRI90 must be changed
C to the one where the ccirnn.asc and ursinn.asc files are.
C
      DO IJF=1,12
        JF(IJF) = .TRUE.
      END DO

      JF(5) = .FALSE.
      JF(12) = .FALSE.

      JMAG = 0
      RZ12 = -F107A
      IDAY = IDATE - IDATE/1000*1000
      MMDD = -IDAY
      CALL IRI90(JF,JMAG,GLAT,GLONG,RZ12,MMDD,STL,Z,JMAX,
     >           'iri/',OUTF,OARR)
      DO J=1,JMAX
        ZE(J) = OUTF(1,J) / 1.E6
        IF (ZE(J) .LT. 100.) ZE(J) = 100.
        ZTI(J) = OUTF(3,J)
        IF (ZTI(J) .LT. ZTN(J)) ZTI(J) = ZTN(J)
        ZTE(J) = OUTF(4,J)
        IF (ZTE(J) .LT. ZTN(J)) ZTE(J) = ZTN(J)
        ZXDEN(3,J) = ZE(J) * OUTF(5,J)/100.
        ZXDEN(6,J) = ZE(J) * OUTF(8,J)/100.
       ZXDEN(7,J) = ZE(J) * OUTF(9,J)/100.
      END DO
C
C
C Fill altitude array and initialize N(2D):
C
      DO J=1,JMAX
        ZZ(J)   = Z(J) * 1.E5
        ZND(J)  = 0.
      END DO
C
C
C Call GLOW to calculate ionized and excited species, airglow emission
C rates, and vertical column brightnesses:
C
      CALL GLOW
C
C
C Set electron densities to calculated values below 200 km, constant
C above:
C
      DO J=JMAX,1,-1
        IF (Z(J) .GT. 200.01) J200=J-1
      END DO
C
      DO J=1,JMAX
        IF (J .LE. J200) ZE(J)=ECALC(J)
        IF (J .GT. J200) ZE(J)=ECALC(J200)
        IF (ZE(J) .LT. 100.) ZE(J) = 100.
      END DO
C
C
C Call GLOW again:
C
      CALL GLOW
C
C
C Output section:
C
      SZAD = SZA * 180. / PI
      DIPD = DIP * 180. / PI

C      write (6,444) IDATE, UT, GLAT, GLONG, F107, F107A, AP
C  444 FORMAT (' Date=',i5,' UT=',f6.0,' Lat=',f5.1,' Lon=',f6.1,
C     >        ' F107=',f4.0,' F107A=',f4.0,' Ap=',f4.0)
C      WRITE (6,445) SZAD, STL, DIPD, EFRAC, IERR
C  445 FORMAT (' SZA=',F5.1,' LST=',F5.2,' Dip=',F5.1,
C     >        ' Ec=',F6.3,' Ie=',I1)
C
C Output photoionization, electron impact ionization,
C electron density, and ion densities:
C
C     write (6,690)
C  690 format ('   Z    Photoion   EIion    Ecalc     O+(2P)    ',
C     >        'O+(2D)    O+(4S)     N+         N2+       O2+       NO+')
C    >        '     O        O2         N2        NO')
      do j=1,jmax
        do i=1,nmaj
          tpi(i) = 0.
          do ns=1,nst
            tpi(i) = tpi(i) + photoi(ns,i,j)
          End Do
        End Do
        totpi = tpi(1) + tpi(2) + tpi(3) + phono(1,j)
        totsi = sion(1,j) + sion(2,j) + sion(3,j)

        Pypi(j) = totpi
        Pysi(j) = totsi
C       write (6,730) z(j),totpi,totsi,ecalc(j),(zxden(i,j),i=1,7)
C    >                zo(j),zo2(j),zn2(j),zno(j)
C  730   format (1x, 0p, f5.1, 1p, 14e10.2)
      End Do

      Pyecalc = ecalc

      Pyion(:,1:7) = transpose(zxden(1:7,:))
      Pyion(:,8)=zo; Pyion(:,9)=zo2; Pyion(:,10)=zn2; Pyion(:,11)=zno

      Pyisr(:,1)=ZE; Pyisr(:,2)=ZTE; Pyisr(:,3)=ZTI
C
C
C Output selected volume emission rates and column brightnesses:
C
C      write (6,780)
C  780 format ('   z     3371   4278   5200   5577   6300',
C     >        '   7320  10400   3466   7774   8446')
C      write (6,790) (z(j), (zeta(iw,j),iw=1,10), j=1,jmax)
C  790 format (1x, f5.1, 10f7.1)
C      write (6,795)  (vcb(iw),iw=1,10)
C  795 format (' VCB:',11f7.0)
C
C
C     CALL ROUT('rt.out',EF,EZ,ITAIL,FRACO,FRACO2,FRACN2)
C
      END SUBROUTINE AURORA
