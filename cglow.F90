module cglow

! Defines shared variables for the GLOW model.
! Replaces the function of common /CGLOW/ and header file glow.h in older
! versions of the model.

! Ben Foster, 1/15.

  implicit none
  save

! Array dimensions:
! JMAX    number of altitude levels
! NBINS   number of energetic electron energy bins
! LMAX    number of wavelength intervals for solar flux
! NMAJ    number of major species
! NEX     number of ionized/excited species
! NW      number of airglow emission wavelengths
! NC      number of component production terms for each emission
! NST     number of states produced by photoionization/dissociation
! NEI     number of states produced by electron impact
! NF      number of types of auroral fluxes

! This replaces glow.h:
!
! JMAX (number of vertical levels) is set to 64 for tgcm runs
! If not a tgcm run, JMAX is set in the namelist input file
!
  integer :: JMAX
!
! These don't change:
!
  integer,parameter :: NBINS=190
  integer,parameter :: LMAX=123
  integer,parameter :: NMAJ=3
  integer,parameter :: NEI=10
  integer,parameter :: NEX=20
  integer,parameter :: NW=20
  integer,parameter :: NC=10
  integer,parameter :: NST=6
  integer,parameter :: NF=4

!  
! This replaces common block /CGLOW/:
!
!     COMMON /CGLOW/
!    >    IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM,
!    >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
!    >    ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
!    >    ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
!    >    ZTN(JMAX), ZTI(JMAX), ZTE(JMAX),
!    >    PHITOP(NBINS), EFLUX(NF), EZERO(NF),
!    >    SZA, DIP, EFRAC, IERR,
!    >    ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX),
!    >    WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
!    >    ENER(NBINS), DEL(NBINS),
!    >    PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
!    >    PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX),
!    >    QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX),
!    >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
!    >    EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX),
!    >    ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW)

  integer :: idate,iscale,jlocal,kchem,ierr
  real    :: ut,glat,glong,f107,f107a,f107p,ap,ef,ec
  real    :: HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC, SZA, DIP, EFRAC

  real,dimension(nf)            :: EFLUX,EZERO
  real,dimension(nw)            :: VCB

  real,allocatable,dimension(:) ::             &                   ! (jmax)
    ZZ, ZO, ZN2, ZO2, ZNO, ZNS, ZND, ZRHO, ZE, &
    ZTN, ZTI, ZTE, QTI, EHEAT, TEZ, ECALC

  real,allocatable,dimension(:)     :: PHITOP, ENER, DEL           ! (nbins) 
  real,allocatable,dimension(:)     :: WAVE1, WAVE2, SFLUX         ! (lmax)
  real,allocatable,dimension(:,:)   :: PESPEC, SESPEC, UFLX, DFLX  ! (nbins,jmax)
  real,allocatable,dimension(:,:)   :: ZMAJ, ZCOL, AURI, PIA, SION ! (nmaj,jmax)
  real,allocatable,dimension(:,:)   :: PHONO                       ! (nst,jmax)
  real,allocatable,dimension(:,:)   :: ZXDEN                       ! (nex,jmax)
  real,allocatable,dimension(:,:)   :: ZETA                        ! (nw,jmax)
  real,allocatable,dimension(:,:,:) :: PHOTOI, PHOTOD              ! (nst,nmaj,jmax)
  real,allocatable,dimension(:,:,:) :: AGLW                        ! (nei,nmaj,jmax)
  real,allocatable,dimension(:,:,:) :: ZCETA                       ! (nc,nw,jmax)

  contains

!-----------------------------------------------------------------------

  subroutine cglow_init

    allocate        &
      (ZZ   (jmax), &
       ZO   (jmax), &
       ZN2  (jmax), &
       ZO2  (jmax), &
       ZNO  (jmax), &
       ZNS  (jmax), &
       ZND  (jmax), &
       ZRHO (jmax), &
       ZE   (jmax), &
       ZTN  (jmax), &
       ZTI  (jmax), &
       ZTE  (jmax), &
       QTI  (jmax), &
       EHEAT(jmax), &
       TEZ  (jmax), &
       ECALC(jmax))

    allocate            &
      (ZXDEN(nex,jmax), &
       ZETA(nw,jmax),   &
       ZCETA(nc,nw,jmax))

    allocate          &
      (PHITOP(nbins), &
       ENER  (nbins), &
       DEL   (nbins))

    allocate          &
      (WAVE1(lmax),   &
       WAVE2(lmax),   &
       SFLUX(lmax))

    allocate               &
      (PESPEC(nbins,jmax), &
       SESPEC(nbins,jmax), &
       UFLX  (nbins,jmax), &
       DFLX  (nbins,jmax))

!    >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),

    allocate &
      (ZMAJ(nmaj,jmax), &
       ZCOL(nmaj,jmax), &
       AURI(nmaj,jmax), &
       PIA (nmaj,jmax), &
       SION(nmaj,jmax)) 

    allocate &
      (AGLW  (nei,nmaj,jmax), &
       PHOTOI(nst,nmaj,jmax), &
       PHOTOD(nst,nmaj,jmax), &
       PHONO(nst,jmax))
    
  end subroutine cglow_init
!-----------------------------------------------------------------------
end module cglow
