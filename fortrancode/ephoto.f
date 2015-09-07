C Subroutine EPHOTO
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Adapted from Banks & Nagy 2-stream input code by Stan Solomon, 6/88
C Modified to handle Auger electrons, Stan Solomon, 7/90
C Reads cross sectons from files (for 1-nm bins), Scott Bailey, ~1994
C Modified bin structure, fixed CIII problem, Stan Solomon, 12/2000
C Corrected additional Auger problem, Liying Qian, 11/2002
C Converged above three branches, Stan Solomon, 3/2005
C Removed LIMIN, wavelength loop now runs from 1 to LMAX, SCS, 3/2005
C
C This subroutine calculates photoionization, rates, certain
C photodissociative excitation rates, and the photoelectron production
C spectrum as a function of altitude.  Uses continuously variable energy
C grid.  3 major species: O, O2, N2; NO is treated as a minor (non-
C absorbing) specie.
C
C Supplied by calling routine:
C WAVE1   wavelength array, upper bound; Angstroms
C WAVE2   wavelength array, lower bound; Angstroms
C SFLUX   solar flux array; photons cm-2 sec-1
C ZZ      altitude array; cm above earth
C ZMAJ    density array for species O, O2, N2, altitude; cm-3
C ZNO     density of NO at each altitude; cm-3
C ZCOL    slant column density for species O, O2, N2, altitude; cm-2
C ENER    energy grid for photoelectrons; eV
C DEL     array of energy grid increments; eV
C
C Calculated by subroutine:
C PESPEC  photoelectron production spectrum for each altitude; cm-3 s-1
C PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
C PHOTOD  photodissoc./exc. rates for state, species, alt.; cm-3 s-1
C PHONO   photoionization/dissoc./exc. rates for NO; cm-3 s-1
C
C Other definitions:
C DSPECT  ionization rate in particular wavelength bin; cm-3 s-1
C TAU     optical depth, dimensionless
C FLUX    solar flux at altitude; cm-2 s-1
C SIGABS  photoabsorption cross sections, O, O2, N2; cm2
C SIGION  photoionization cross sections, O, O2, N2; cm2
C SIGAO, SIGAO2, SIGAN2, SIGIO, SIGIO2, SIGIN2; cross sect. data arrays
C NNN     number of states for each species
C TPOT    ionization potentials for each species, state; eV
C PROB    branching ratios for each state, species, and wavelength bin:
C         O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
C         O2+ states: X, a+A, b, dissoc.
C         N2+ states: X, A, B, C, F, dissoc.
C PROBO, PROBO2, PROBN2; branching ratio data arrays
C BSO2    yield of O(1S) from dissociation of O2
C EPSIL1  energy loss lower bound for state, species, wavelength; eV
C EPSIL2  energy loss upper bound for state, species, wavelength; eV
C SIGNO   NO photoionization xsect at Ly-alpha
C AUGE    Mean energy of Auger electrons for each species; eV
C AUGL    Wavelength threshold for Auger electrons; Angstroms
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
      SUBROUTINE EPHOTO
C
C
      INCLUDE 'glow.h'
      PARAMETER (NMAJ=3)
      PARAMETER (NEX=20)
      PARAMETER (NW=20)
      PARAMETER (NC=10)
      PARAMETER (NST=6)
      PARAMETER (NEI=10)
      PARAMETER (NF=4)
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
      DIMENSION DSPECT(JMAX), FLUX(LMAX,JMAX), NNN(NMAJ),
     >          SIGION(NMAJ,LMAX), SIGABS(NMAJ,LMAX),
     >          TPOT(NST,NMAJ), PROB(NST,NMAJ,LMAX),
     >          EPSIL1(NST,NMAJ,LMAX), EPSIL2(NST,NMAJ,LMAX),
     >          SIGAO(LMAX), SIGAO2(LMAX), SIGAN2(LMAX),
     >          SIGIO(LMAX), SIGIO2(LMAX), SIGIN2(LMAX),
     >          PROBO(NST,LMAX), PROBO2(NST,LMAX), PROBN2(NST,LMAX),
     >          BSO2(LMAX), AUGE(NMAJ), AUGL(NMAJ), TAU(LMAX),
     >          RION(LMAX,NMAJ,JMAX)
C
      SAVE SIGION, SIGABS, PROB, EPSIL1, EPSIL2
C
      DATA SIGNO/2.0 E-18/, NNN/5,4,6/, IFIRST/1/
C
      DATA TPOT/13.61, 16.93, 18.63, 28.50, 40.00,  0.00,
     >          12.07, 16.10, 18.20, 20.00,  0.00,  0.00,
     >          15.60, 16.70, 18.80, 30.00, 34.80, 25.00/
C
      DATA BSO2/12*0.,.01,.03,7*.10,8*.07,5*.03,5*.01,84*0./
C
      DATA AUGE/500., 500., 360./, AUGL/24., 24., 33./
C
C
C NB - absorption and ionization cross sections are multiplied by 1.E-18
C on first call.
C
C
C First time only: Read cross section data from files, convert to cm2,
C calculate energy losses:
C
      IF (IFIRST .EQ. 1) THEN
      IFIRST = 0
C
      open(unit=1,file='ephoto_xn2.dat',status='old')
          read(1,*)
          read(1,*)
          read(1,*)
          read(1,*)
          Do l=lmax,1,-1
              read(1,*) aa,bb,(probn2(n,l),n=1,nst),sigin2(l),sigan2(l)
              bso2(l)=0.0
          End Do
      close(1)
C
      open(unit=1,file='ephoto_xo2.dat',status='old')
          read(1,*)
          read(1,*)
          read(1,*)
          read(1,*)
          Do l=lmax,1,-1
             read(1,*) aa,bb,(probo2(n,l),n=1,nst),sigio2(l),sigao2(l)
          End Do
      close(1)
C 
      open(unit=1,file='ephoto_xo.dat',status='old')
          read(1,*)
          read(1,*)
          read(1,*)
          read(1,*)
          do l=lmax,1,-1
             read(1,*) aa,bb,(probo(n,l),n=1,nst),sigio(l),sigao(l)
          End Do
      close(1)
C
        DO L=1,LMAX
            SIGABS(1,L) = SIGAO(L)  * 1.E-18
            SIGABS(2,L) = SIGAO2(L) * 1.E-18
            SIGABS(3,L) = SIGAN2(L) * 1.E-18
            SIGION(1,L) = SIGIO(L)  * 1.E-18
            SIGION(2,L) = SIGIO2(L) * 1.E-18
            SIGION(3,L) = SIGIN2(L) * 1.E-18
        End Do
C
        DO L=1,LMAX
            DO K=1,NST
                PROB(K,1,L) = PROBO(K,L)
                PROB(K,2,L) = PROBO2(K,L)
                PROB(K,3,L) = PROBN2(K,L)
            End Do
        End Do
C
        DO L=1,LMAX 
            DO I=1,NMAJ 
                DO K=1,NNN(I) 
                    EPSIL1(K,I,L)=12397.7/WAVE1(L)-TPOT(K,I) 
                    EPSIL2(K,I,L)=12397.7/WAVE2(L)-TPOT(K,I) 
                    IF (WAVE1(L) .LE. AUGL(I)) THEN
                      EPSIL1(K,I,L) = EPSIL1(K,I,L) - AUGE(I)
                      EPSIL2(K,I,L) = EPSIL2(K,I,L) - AUGE(I)
                    ENDIF
                End Do
            End Do
        End Do 
C
      ENDIF
C
C
C Zero arrays:
C
      DO J=1,JMAX
          DO I=1,NMAJ
              DO K=1,NST
                  PHONO(K,J) = 0.
                  PHOTOI(K,I,J) = 0.
                  PHOTOD(K,I,J) = 0.
              End Do
          End Do
          DO M=1,NBINS
            PESPEC(M,J) = 0.
          End Do
      End Do
C
C
C Calculate attenuated solar flux at all altitudes and wavelengths:
C
      DO 200 L=1,LMAX 
      DO 200 J=1,JMAX
      TAU(L)=0. 
      DO I=1,NMAJ 
      	TAU(L)=TAU(L)+SIGABS(I,L)*ZCOL(I,J) 
      End Do
      IF (TAU(L) .LT. 20.) THEN
        FLUX(L,J)=SFLUX(L)*EXP(-TAU(L)) 
      ELSE
        FLUX(L,J) = 0.0
      ENDIF
C
C
C Calculate SRC photodissociation of O2, dissociative excitation of
C O(1S), photodissociation of N2, and photoionization of NO by solar
C Ly-alpha:
C
      IF (WAVE1(L) .LT. 1751. AND. WAVE2(L) .GT. 1349.)
     >   PHOTOD(1,2,J) = PHOTOD(1,2,J)+ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J)
      PHOTOD(2,2,J) = PHOTOD(2,2,J) + ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J)
     >                                * BSO2(L)
      PHOTOD(1,3,J) = PHOTOD(1,3,J) +
     >                ZMAJ(3,J)*(SIGABS(3,L)-SIGION(3,L))*FLUX(L,J)
      IF (WAVE1(L) .LT. 1221. .AND. WAVE2(L) .GT. 1209.)
     >   PHONO(1,J) = PHONO(1,J) + ZNO(J)*SIGNO*FLUX(L,J)
  200 CONTINUE
C
C
C
C Calculate ionization rates and photoelectron production:
C
C
C Loop over wavelengths:
C
      DO 400 L=1,LMAX
C
C
C Loop over species:
C
      DO 350 I=1,NMAJ 
C
C
C Calculate total ionization rates for all species and altitudes:
C
      DO J=1,JMAX
        RION(L,I,J)=ZMAJ(I,J)*SIGION(I,L)*FLUX(L,J)
      End Do
C
C
C Loop over states:
C
      DO 300 K=1,NNN(I) 
C
      E1= EPSIL1(K,I,L) 
      E2= EPSIL2(K,I,L) 
      IF (E2 .LT. 0.) GO TO 300 
      IF (E1 .LT. 0.) E1=0. 
C
C
C Calculate state-specific ionization rates at all altitudes:
C
      DO J=1,JMAX
        DSPECT(J) = RION(L,I,J)*PROB(K,I,L) 
        PHOTOI(K,I,J) = PHOTOI(K,I,J) + DSPECT(J)
      End Do
C
C
C Find box numbers M1, M2 corresponding to energies E1, E2:
C
      CALL BOXNUM (E1, E2, M1, M2, R1, R2, NBINS, DEL, ENER) 
      IF (M1 .GT. NBINS) GOTO 300
C
C
C Fill the boxes from M1 to M2 at all altitudes:
C 
      Y = E2 - E1 
      DO N=M1,M2
          IF (M1 .EQ. M2) THEN
            FAC = 1.
          ELSE
                IF (N .EQ. M1) THEN
                  FAC = (R1-E1) / Y
                ELSE
                  IF (N .EQ. M2) THEN
                    FAC = (E2-R2) / Y
                  ELSE
                    FAC = DEL(N) / Y
                  ENDIF
                ENDIF
          ENDIF
          DO J=1,JMAX
              PESPEC(N,J) = PESPEC(N,J) + DSPECT(J) * FAC
          End Do
      End Do
C
  300 CONTINUE 
C
C
C Generate Auger electrons if energy is sufficient:
C
      IF (WAVE1(L) .LE. AUGL(I)) THEN
        E1 = AUGE(I)
        E2 = AUGE(I)
        CALL BOXNUM (E1, E2, M1, M2, R1, R2, NBINS, DEL, ENER) 
        IF (M1.GT.NBINS .OR. M2.GT.NBINS) GOTO 350
        DO J=1,JMAX
          PESPEC(M1,J) = PESPEC(M1,J) + RION(L,I,J)
        End Do
      ENDIF
C     
C
  350 CONTINUE
C
  400 CONTINUE

      END SUBROUTINE EPHOTO
C
C
C
C
      SUBROUTINE BOXNUM (E1, E2, M1, M2, R1, R2, NBINS, DEL, ENER)
      implicit none
C This subroutine finds the box numbers corresponding to
C energies E1 and E2, and calls them M1 and M2
C 
C R1 is the upper edge of the lower box, R2 is the lower edge of the
C upper box.
C
      Integer,Intent(in) :: NBINS
      Real,Intent(in)    :: DEL(NBINS), ENER(NBINS)
      Real,Intent(out)   :: E1,E2,R1,R2
      Integer,Intent(out):: M1,M2
      
      Integer I
C
      DO I=1,NBINS
          IF (E1 .LT. ENER(I)+DEL(I)/2.) GOTO 200
      End Do
      M1 = NBINS+1
      RETURN
C
  200 M1 = I
      R1 = ENER(I) + DEL(I)/2.
C
      DO I=1,NBINS
      	IF (E2 .LT. ENER(I)+DEL(I)/2.) GOTO 400
      End Do
      M2 = NBINS
      R2 = E2 - DEL(NBINS)
      RETURN
C
  400 M2 = I
      R2 = ENER(I) - DEL(I)/2.
C
C
      END Subroutine BOXNUM
