C Subroutine ETRANS
C
! M.H.: this appears to calculate electron precipitation impact outcomes on various species
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Banks & Nagy 2-stream electron transport code
C Adapted by Stan Solomon, 1986, 1988
C Uses variable altitude and energy grids
C LMAX obtained from glow.h, SMB, 1994
C Updated comments and removed artifacts, SCS, 2005
C
C Subroutine EXSECT called first time only to calculate electron impact
C cross sections.
C
C Definitions:
C COMMON /CGLOW/:  see subroutine GLOW
C COMMON /CXSECT/: see subroutine EXSECT
C COMMON /CXPARS/: see subroutine EXSECT
C PSI    first term of parabolic d.e., = 1
C ALPHA  second term "; cm-1
C BETA   third term  "; cm-2
C GAMA  forth term  "; cm-4 s-1 eV-1
C DELZ   altitude increments; cm
C DEL2   sum of altitude increment and next higher increment; cm
C DELA   average of "
C DELP   product of DELA and next higher DELZ
C DELM   product of DELA and DELZ
C DELS   product of DELZ and next higer DELZ
C DEN    dummy array for transfer of calculated downward flux
C FAC    factor for extrapolating production rate, = 0
C PROD   sum of photo- and secondary electron production; cm-3 s-1 eV-1
C EPROD  energy of "; eV cm-3
C T1     elastic collision term; cm-1
C T2     elastic + inelastic collision term; cm-1
C TSA    total energy loss cross section for each species; cm2
C PRODUP upward   cascade + secondary production; cm-3 s-1 eV-1
C PRODWN downward cascade + secondary production; cm-3 s-1 eV-1
C PHIUP  upward   flux; cm-2 s-1 eV-1
C PHIDWN downward flux; cm-2 s-1 eV-1
C TSIGNE thermal electron collision term; cm-1
C SECION total ionization rate; cm-3 s-1
C SECP   secondary electron production; cm-3 s-1 eV-1
C R1     ratio term for calculating upward flux; cm-2 s-1 eV-1
C EXPT2  exponential term for calculating upward flux
C PRODUA collection array for calculating PRODUP; cm-3 s-1 eV-1
C PRODDA  "                               PRODWN
C PHIINF downward flux at top of atmos., divided by AVMU; cm-2 s-1 eV-1
C POTION ionizaition potential for each species; eV
C AVMU   cosine of the average pitch angle
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
C NF      number of available types of auroral fluxes
C
C
      SUBROUTINE ETRANS
!      use cglow, only: nmaj,jmax,nw,nst,nbins,nei,nf,lmax,nc,nex
! *************************
! TODO using "implicit none" causes this file to give error with f2py/f2py3
!constructing wrapper function "aurora"...
!         pyion,pyecalc,pypi,pysi,pyisr = aurora(z,pyidate,pyut,pyglat,pyglong,pyf107a,pyf107,pyf107p,pyap,pyphitop)
!{}
!analyzevars: charselector={'len': '4'} unhandled.analyzevars: charselector={'len': '4'} unhandled.analyzevars: charselector={'len': '4'} unhandled.getctype: No C-type found in "{}", assuming void.
!
! issue with common block? Leave it alone till Stan has new F90 code available.
! ************************

      use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, 
     &                                          stderr=>error_unit
!      implicit none
      INCLUDE 'cglow.h'

      logical isfinite

      integer  IDATE, ISCALE, JLOCAL, KCHEM, IERR, IIMAXX(NBINS)

      real   UT, GLAT, GLONG,
     >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
     >    ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
     >    ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
     >    ZTN(JMAX), ZTI(JMAX), ZTE(JMAX),
     >    PHITOP(NBINS), EFLUX(NF), EZERO(NF),
     >    SZA, DIP, EFRAC,
     >    ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX),
     >    WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
     >    ENER(NBINS), DEL(NBINS),
     >    PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
     >    PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX),
     >    QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX),
     >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
     >    EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX),
     >    ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW),
     &    SIGS(NMAJ,NBINS), PE(NMAJ,NBINS), PIN(NMAJ,NBINS),
     >                SIGA(NMAJ,NBINS,NBINS), SEC(NMAJ,NBINS,NBINS),
     >                SIGEX(NEI,NMAJ,NBINS), SIGIX(NEI,NMAJ,NBINS),
     &    WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ),
     >                ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ),
     >                THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ),
     >                TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ),
     >                GAMS(NEI,NMAJ), GAMB(NEI,NMAJ),
     &    ALPHA(JMAX), BETA(JMAX), GAMA(JMAX), PSI(JMAX),
     >                DELZ(JMAX), DEL2(JMAX), DELA(JMAX), DELP(JMAX),
     >                DELM(JMAX), DELS(JMAX), DEN(JMAX), FAC

      real    PROD(JMAX), EPROD(JMAX), T1(JMAX), T2(JMAX), TSA(NMAJ),
     >          PRODUP(JMAX,NBINS), PRODWN(JMAX,NBINS),
     >          PHIUP(JMAX), PHIDWN(JMAX), TSIGNE(JMAX), TAUE(JMAX),
     >          SECION(JMAX), SECP(NMAJ,JMAX), R1(JMAX), EXPT2(JMAX),
     >          PRODUA(JMAX), PRODDA(JMAX), PHIINF(NBINS)

      real  APROD,DAG,EDEP,EET,ein,eout,epe,ephi,et,fluxj,phiout,
     &     rmusin, sindip
      integer  i,ib,ibb,ii,im,iq,iv,j,jj,jjj4,k,kk,ll,n

      COMMON /CGLOW/ IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM,
     >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
     >    ZZ, ZO, ZN2, ZO2, ZNO, ZNS, ZND, ZRHO, ZE,
     >    ZTN, ZTI, ZTE, PHITOP, EFLUX, EZERO, SZA, DIP, EFRAC, IERR,
     >    ZMAJ, ZCOL, WAVE1, WAVE2, SFLUX, ENER, DEL, PESPEC, SESPEC,
     >    PHOTOI, PHOTOD, PHONO, QTI, AURI, PIA, SION,
     >    UFLX, DFLX, AGLW, EHEAT, TEZ, ECALC, ZXDEN, ZETA, ZCETA, VCB


      COMMON /CXSECT/ SIGS, PE, PIN, SIGA, SEC, SIGEX, SIGIX, IIMAXX

      COMMON /CXPARS/ WW, AO, OMEG, ANU, BB, AUTO,THI, AK, AJ,
     >                TS, TA, TB, GAMS, GAMB

      COMMON /CIMPIT/ ALPHA, BETA, GAMA, PSI,DELZ, DEL2, DELA, DELP,
     >                DELM, DELS, DEN, FAC

      integer :: IFIRST=1
      real,parameter :: AVMU=0.5

      real potion(nmaj)
      DATA POTION/16.,16.,18./

      IERR = 0
      FAC = 0.
      SINDIP = SIN(DIP)
      RMUSIN = 1. / SINDIP / AVMU

C First call only:  calculate cross-sectons:

      IF (IFIRST .EQ. 1) THEN
        CALL EXSECT (ENER, DEL)
        IFIRST = 0
      ENDIF
C
C
C Zero variables:
C
      DO 100 II=1,JMAX
        DO 100 IB=1,NMAJ
          DO 100 IBB=1,NEI
            AGLW(IBB,IB,II) = 0.0
  100 CONTINUE
C
      PSI(1)   = 1.
      ALPHA(1) = 0.
      BETA(1) = 0.
      GAMA(1) = 0.
      PHIOUT = 0.0
C
      DO I = 1, JMAX
        EHEAT(I) = 0.0
        EPROD(I) = 0.0
        SECION(I) = 0.0
        DO N = 1, NMAJ
          SION(N,I) = 0.0
        End Do
      End Do
C
      DO 500 JJ = 1, NBINS
        DO 500 I = 1, JMAX
          PRODUP(I,JJ) = 1.0E-20
          PRODWN(I,JJ) = 1.0E-20
  500 CONTINUE
C
C
C Divide downward flux at top of atmos. by average pitch angle cosine:
C
      DO  J=1,NBINS
        PHIINF(J) = PHITOP(J) / AVMU
        if(debug .and. .not.isfinite(phiinf(j)))
     &             error stop 'etrans: nonfinite PhiInf'
      End Do
C
C
C Calcualte delta z's:
C
      DELZ(1) = ZZ(2)-ZZ(1)
      DO I=2,JMAX
        DELZ(I) = ZZ(I)-ZZ(I-1)
      End Do

      DO I=1,JMAX-1
        DEL2(I) = DELZ(I)+DELZ(I+1)
        DELA(I) = DEL2(I)/2.
        DELP(I) = DELA(I)*DELZ(I+1)
        DELM(I) = DELA(I)*DELZ(I)
        DELS(I) = DELZ(I)*DELZ(I+1)
      End Do

      DEL2(JMAX) = DEL2(JMAX-1)
      DELA(JMAX) = DELA(JMAX-1)
      DELP(JMAX) = DELP(JMAX-1)
      DELM(JMAX) = DELP(JMAX-1)
      DELS(JMAX) = DELS(JMAX-1)
C
C
C
C Top of Energy loop:
C
      DO 1100 J=NBINS,1,-1
C
C
C Calculate production:
C
      DO I = 1, JMAX
          PROD(I) = (PESPEC(J,I)+SESPEC(J,I)) * RMUSIN / DEL(J)
          EPROD(I) = EPROD(I) + PROD(I) * ENER(J) * DEL(J) / RMUSIN
      End Do
C
C
C Total energy loss cross section for each species:
C
      DO 740 I = 1, NMAJ
        TSA(I) = 0.0
  740 CONTINUE
      IF (J .GT. 1) THEN
        DO 760 K = 1, J-1
          DO 760 I = 1, NMAJ
            TSA(I) = TSA(I) + SIGA(I,K,J) * (DEL(J-K)/DEL(J))
  760   CONTINUE
      ELSE
        DO I=1,NMAJ
          TSA(I) = TSA(I) + SIGA(I,1,J) + 1.E-18
        End Do
      ENDIF
C
C
C Thermal electron energy loss:
C
      JJJ4 = J - 1
      IF (J .EQ. 1) JJJ4 = 1
      DAG = ENER(J) - ENER(JJJ4)
      IF (DAG .LE. 0.0) DAG = DEL(1)
C
      DO 800 I = 1, JMAX
        ET = 8.618E-5 * ZTE(I)
        EET = ENER(J) - ET
        IF (EET .LE. 0.0) GO TO 780
        TSIGNE(I) = ((3.37E-12*ZE(I)**0.97)/(ENER(J)**0.94))
     >              * ((EET)/(ENER(J) - (0.53*ET))) ** 2.36
        GO TO 790
  780   TSIGNE(I) = 0.0
  790   CONTINUE
        TSIGNE(I) = TSIGNE(I) * RMUSIN / DAG
  800 CONTINUE
C
C
C Collision terms:
C
      DO 870 I = 1, JMAX
        T1(I) = 0.0
        T2(I) = 0.0
        DO IV = 1, NMAJ
          T1(I) = T1(I) + ZMAJ(IV,I) * SIGS(IV,J) * PE(IV,J)
          T2(I) = T2(I) + ZMAJ(IV,I) * (SIGS(IV,J)*PE(IV,J) + TSA(IV))
        End Do
        T1(I) = T1(I) * RMUSIN
        T2(I) = T2(I) * RMUSIN + TSIGNE(I)
  870 CONTINUE
C
C
C Bypass next section if local calculation was specified:
C
      IF (JLOCAL .EQ. 1) GO TO 910
C
C
C Solve parabolic d.e. by Crank-Nicholson method to find downward flux:
C
      DO I = 2, JMAX-1
        PSI(I) = 1.
        ALPHA(I) = (T1(I-1) - T1(I+1)) / (DEL2(I) * T1(I))
        BETA(I) = T2(I) * (T1(I+1) - T1(I-1)) / (T1(I) * DEL2(I))
     >            - (T2(I+1) - T2(I-1)) / DEL2(I)
     >            - T2(I)**2 + T1(I)**2
        IF (PROD(I) .LT. 1.E-30) PROD(I) = 1.E-30
        IF (PRODWN(I,J) .LT. 1.E-30) PRODWN(I,J) = 1.E-30
        GAMA(I) = (PROD(I)/2.0)
     >             * (-T1(I) - T2(I) - ALPHA(I)
     >                - (PROD(I+1) - PROD(I-1))/PROD(I)/DEL2(I))
     >             + PRODWN(I,J)
     >             * (-ALPHA(I) - T2(I)
     >                - (PRODWN(I+1,J) - PRODWN(I-1,J))
     >                   /PRODWN(I,J)/DEL2(I))
     >             - PRODUP(I,J) * T1(I)
        !if (.not.isfinite(GAMA(I))) GAMA(I) = 0. 
        if (debug .and. .not.isfinite(GAMA(I))) then
         write (stderr,*),'etrans.f: GAMA PRODWNn1 PRODWN PRODWNp1',
     &     ' PRODUP PRODn1 PROD PRODp1 T1 T2 ALPHA DEL2'
         write(stderr,*), GAMA(I),PRODWN(I-1,J),PRODWN(I,J),
     &     PRODWN(I+1,J), PRODUP(I,J),PROD(I-1),PROD(I),PROD(I+1),
     &     T1(I),T2(I),ALPHA(I), DEL2(I)
         error stop 'etran: non-finite GAMA'
        end if
      End DO

      IF (ABS(BETA(2)) .LT. 1.E-20) THEN
        BETA(2) = 1.E-20
        IERR = 2
      ENDIF
      PHIDWN(2) = GAMA(2) / BETA(2)
      DEN(1) = PHIDWN(2)
      FLUXJ = PHIINF(J)
      CALL IMPIT(FLUXJ) !computes DEN via COMMON array
      DO I = 1, JMAX
        PHIDWN(I) = DEN(I)
        if(debug .and. phidwn(i).gt.1e30) then
            write(stderr,*) 'GAMA(i), beta(i)',gama(i),beta(i)
          error stop 'etrans: very large PHIDWN (jlocal != 1)'
        endif
      End Do
C
C
C Apply lower boundary condition: PHIUP=PHIDWN.  Should be nearly zero.
C Then integrate back upward to calculate upward flux:
C
      PHIUP(1) = PHIDWN(1)
      DO 900 I = 2, JMAX
        R1(I) = (T1(I)*PHIDWN(I) + (PROD(I)+2.*PRODUP(I,J))/2.) / T2(I)
        TAUE(I) = T2(I)*DELZ(I)
        IF (TAUE(I) .GT. 60.) TAUE(I)=60.
        EXPT2(I) = EXP(-TAUE(I))
  900 CONTINUE
      DO I=2,JMAX
        PHIUP(I) = R1(I) + (PHIUP(I-1)-R1(I)) * EXPT2(I)        
        if (debug .and. .not.isfinite(phiup(i))) 
     &           error stop 'etrans: nonfinite PHIUP  (jlocal != 1)'
      End DO

      GOTO 930
C
C
C Local calculation only  (if jlocal=1):
C
  910 CONTINUE
      DO I = 1, JMAX
        IF (T2(I) .LE. T1(I)) THEN
          IERR = 1
          T2(I) = T1(I) * 1.0001
        ENDIF
        PHIUP(I) = (PROD(I)/2.0 + PRODUP(I,J)) / (T2(I) - T1(I))
        PHIDWN(I) = (PROD(I)/2.0 + PRODWN(I,J)) / (T2(I) - T1(I))
        if (debug .and. .not. isfinite(phiup(i))) 
     &       error stop 'etrans: nonfinite PHIUP  (jlocal=1)'
        if (debug .and. phidwn(i).gt.1e30) 
     &       error stop 'etrans: very large PHIDWN (jlocal=1)'
      End Do
C
  930 CONTINUE  ! if jlocal != 1
C
C
C Multiply fluxes by average pitch angle cosine and put in arrays,
C and calculate outgoing electron energy flux for conservation check:
C
      DO I=1,JMAX
        UFLX(J,I) = PHIUP(I) * AVMU
        DFLX(J,I) = PHIDWN(I) * AVMU
      End do
C
      PHIOUT = PHIOUT + PHIUP(JMAX) * DEL(J) * ENER(J)
C
C
C Cascade production:
C
      DO 990 K = 1, J-1
        LL = J - K
        DO 972 I=1,JMAX
          PRODUA(I) = 0.0
          PRODDA(I) = 0.0
  972   CONTINUE
        DO 980 N = 1, NMAJ
          DO 975 I=1,JMAX
            PRODUA(I) = PRODUA(I)
     >                  + ZMAJ(N,I) * (SIGA(N,K,J)*PIN(N,J)*PHIDWN(I)
     >                  + (1. - PIN(N,J))*SIGA(N,K,J)*PHIUP(I))
            PRODDA(I) = PRODDA(I)
     >                  + ZMAJ(N,I) * (SIGA(N,K,J)*PIN(N,J)*PHIUP(I)
     >                  + (1. - PIN(N,J))*SIGA(N,K,J)*PHIDWN(I))
  975     CONTINUE
  980   CONTINUE
        DO 985 I=1,JMAX
          PRODUP(I,LL) = PRODUP(I,LL) + PRODUA(I) * RMUSIN
          PRODWN(I,LL) = PRODWN(I,LL) + PRODDA(I) * RMUSIN
  985   CONTINUE
  990 CONTINUE
C
      KK = J - 1
      IF (KK .LE. 0) GO TO 1020
      DO 1010 I = 1, JMAX
        PRODUP(I,KK) = PRODUP(I,KK) + TSIGNE(I) * PHIUP(I) * (DEL(J)
     >                 / DEL(KK))
        PRODWN(I,KK) = PRODWN(I,KK) + TSIGNE(I) * PHIDWN(I) * (DEL(J)
     >                 / DEL(KK))
 1010 CONTINUE
 1020 CONTINUE
C
C
C Electron heating rate:
C
      DAG = DEL(J)
      DO I = 1, JMAX
        EHEAT(I) = EHEAT(I) + TSIGNE(I) * (PHIUP(I)+PHIDWN(I)) * DAG**2
      End Do
C
C
C Electron impact excitation rates:
C
      DO 1040 II = 1, JMAX
        DO 1040 I = 1, NMAJ
          DO 1040 IBB = 1, NEI
            AGLW(IBB,I,II) = AGLW(IBB,I,II) + (PHIUP(II) + PHIDWN(II))
     >                       * SIGEX(IBB,I,J) * DEL(J) * ZMAJ(I,II)
 1040 CONTINUE
C
C
C Calculate production of secondaries into K bin for energy J bin and
C add to production:
C
      DO 1090 K = 1, IIMAXX(J) ! iimaxx set near exsect.f:424
        DO 1080 N = 1, NMAJ
          DO 1070 I = 1, JMAX
            SECP(N,I) = SEC(N,K,J) * ZMAJ(N,I) * (PHIUP(I) + PHIDWN(I))
!            if (isnan(sec(n,k,j))) stop 'etrans: NaN in SEC'
!            if (isnan(zmaj(n,i))) stop 'etrans: NaN in ZMAJ'
!            if (isnan(phiup(i))) stop 'etrans: NaN in PHIUP'
            if (debug .and. phidwn(i).gt.1e30) 
     &            error stop 'etrans.f: very large PHIDWN'
            if (debug .and. .not.isfinite(secp(n,i)))
     &            error stop 'etrans: nonfinite SECP'
            SION(N,I) = SION(N,I) + SECP(N,I) * DEL(K)
            
            if (debug .and. .not.isfinite(sion(n,i)))  then
             write(stderr,*)'etrans: nonfinite impact ioniz SION.',
     &                 '  k,n,i=',k,n,i
             write(stderr,*) 'del(k)=',del(k)
             write(stderr,*) 'secp(n,i)=',secp(n,i)
             write(stderr,*) 'sion(n,i)=',sion(n,i)
            error stop
            endif

            SECION(I) = SECION(I) + SECP(N,I) * DEL(K)
            PRODUP(I,K) = PRODUP(I,K) + (SECP(N,I)*.5*RMUSIN)
            PRODWN(I,K) = PRODWN(I,K) + (SECP(N,I)*.5*RMUSIN)
 1070     CONTINUE
 1080   CONTINUE
 1090 CONTINUE
C
 1100 CONTINUE
C
C Bottom of Energy loop
C

      DO 1250 I = 1, JMAX
        EHEAT(I) = EHEAT(I) / RMUSIN
 1250 CONTINUE
C
C
C Calculate energy deposited as a function of altitude
C and total energy deposition:
C
      EDEP = 0.
      DO 1270 IM=1,JMAX
        TEZ(IM) = EHEAT(IM)
        DO 1260 II=1,NMAJ
          TEZ(IM) = TEZ(IM) + SION(II,IM)*POTION(II)
          DO 1260 IQ=1,NEI
            TEZ(IM) = TEZ(IM) + AGLW(IQ,II,IM)*WW(IQ,II)
 1260   CONTINUE
        EDEP = EDEP + TEZ(IM) * DELA(IM)
 1270 CONTINUE
C
C
C Calculate energy input, output, and fractional conservation:
C
      EPE = 0.0
      EPHI = 0.0
      DO 1440 I = 2, JMAX
        APROD = SQRT(EPROD(I)*EPROD(I - 1))
        EPE = EPE + APROD * DELZ(I)
 1440 CONTINUE
      DO JJ = 1, NBINS
        EPHI = EPHI + PHIINF(JJ) * ENER(JJ) * DEL(JJ) / RMUSIN
      End Do
      EIN = EPHI + EPE
      PHIOUT = PHIOUT / RMUSIN
      EOUT = EDEP + PHIOUT
      EFRAC = (EOUT - EIN) / EIN

      END SUBROUTINE ETRANS
C
C
C
C
C Subroutine IMPIT solves parabolic differential equation by implicit
C Crank-Nicholson method
C
      SUBROUTINE IMPIT(FLUXJ)
      use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, 
     &                                          stderr=>error_unit
!      use cglow,only: jmax
      Implicit None
      include 'cglow.h'
!Args:
      Real, Intent(In) :: FLUXJ
!Local:
      logical isfinite
      Real fac,dem
      real,dimension(jmax) :: alpha, beta, gama, psi, delz, del2, dela,
     > delp,delm,dels,den,
     > K, L, A, B, C, D
      Integer i,i1,jk,kk

      COMMON /CIMPIT/ ALPHA, BETA, GAMA, PSI, DELZ, DEL2, DELA, DELP,
     >                DELM, DELS, DEN, FAC
C
      I1 = JMAX - 1

      DO I = 1, I1
        A(I) = PSI(I) / DELP(I) + ALPHA(I) / DEL2(I)
        B(I) = -2. * PSI(I) / DELS(I) + BETA(I)
        C(I) = PSI(I) / DELM(I) - ALPHA(I) / DEL2(I)
        D(I) = GAMA(I)
       if(debug .and. .not.isfinite(c(i)))
     &              error stop 'etrans:impit nonfinite C(I)'
       if(debug .and. .not.isfinite(d(i)))
     &              error stop 'etrans:impit nonfinite D(I)'
      End Do

      K(2) = (D(2) - C(2)*DEN(1)) / B(2)
      L(2) = A(2) / B(2)
      if (debug .and. .not.isfinite(k(2))) then
        write(stderr,*) '***********    **********'
        write(stderr,*) 'K(2) D(2) C(2) DEN(1) B(2)',K(2),D(2),C(2),
     &   DEN(1),B(2)
        error stop 'etrans:impit non-finite K(2)'
      end if
!      if (isnan(k(2))) stop 'etrans:impit NaN in K(2)'

      DO I = 3, I1
        DEM = B(I) - C(I) * L(I-1)
        if (debug .and. .not.isfinite(dem)) 
     &        error stop 'etrans:impit nonfinite DEM'
        K(I) = (D(I) - C(I)*K(I-1)) / DEM
        L(I) = A(I) / DEM

        if (debug .and. .not. isfinite(K(i))) then
         write(stderr,*) k
         error stop 'etrans:impit NaN in K(i)'
        end if
!        if (isnan(L(i))) stop'etrans:impit NaN in L(i)'
      End DO   

      DEN(I1) = (K(I1) - L(I1)*FLUXJ) / (1. + L(I1)*FAC)
!      if (isnan(K(i1))) stop'etrans:impit NaN in K(i1)'
!      if (isnan(L(i1))) stop'etrans:impit NaN in L(i1)'
      if(debug .and. .not.isfinite(den(i1)))
     &          error stop 'impit nonfinite DEN(I1)'
      DEN(JMAX) = DEN(I1)
      
      Do KK = 1, JMAX-3
        JK = I1 - KK
        DEN(JK) = K(JK) - L(JK) * DEN(JK + 1)
        if (debug .and. .not.isfinite(den(jk)))
     &    error stop 'etrans:impit: non-finite DEN'

      End Do

      END Subroutine IMPIT


      logical Function isfinite(x)
      implicit none

      real,intent(in) :: x

      if (abs(x) <= huge(x)) then
        isfinite = .true.
      else
        isfinite = .false.
      end if

      end function isfinite
