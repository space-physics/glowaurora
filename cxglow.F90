module cxglow
!
! This module replaces commons CXSECT, CXPARS, CIMPIT.
!
  use cglow,only: JMAX,NBINS,NMAJ,NEI
  implicit none
  save
!--------------------------------------------------------------------------------
!     COMMON /CXSECT/ SIGS(NMAJ,NBINS), PE(NMAJ,NBINS), PIN(NMAJ,NBINS),
!    >                SIGA(NMAJ,NBINS,NBINS), SEC(NMAJ,NBINS,NBINS),
!    >                SIGEX(NEI,NMAJ,NBINS), SIGIX(NEI,NMAJ,NBINS),
!    >                IIMAXX(NBINS)

  real,allocatable,dimension(:,:)   :: sigs,pe,pin ! (nmaj,nbins)
  real,allocatable,dimension(:,:,:) :: sigex,sigix ! (nei,nmaj,nbins)
  real,allocatable,dimension(:,:,:) :: siga,sec    ! (nei,nbins,nbins)

  integer,allocatable,dimension(:)  :: iimaxx      ! (nbins)
!--------------------------------------------------------------------------------
!     COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ),
!    >                ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ),
!    >                THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ),
!    >                TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ),
!    >                GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)

  real,allocatable,dimension(:,:) :: & ! (nei,nmaj)
    ww,ao,omeg,anu,bb,auto,thi,ak,aj,ts,ta,tb,gams,gamb
!--------------------------------------------------------------------------------
!     COMMON /CIMPIT/ ALPHA(JMAX), BETA(JMAX), GAMMA(JMAX), PSI(JMAX),
!    >                DELZ(JMAX), DEL2(JMAX), DELA(JMAX), DELP(JMAX),
!    >                DELM(JMAX), DELS(JMAX), DEN(JMAX), FAC

  real,allocatable,dimension(:) :: & ! (jmax)
    alpha,beta,gamma,psi,delz,del2,dela,delp,delm,dels,den
  real :: fac
!--------------------------------------------------------------------------------
  contains
!--------------------------------------------------------------------------------
  subroutine cxglow_init

    allocate             &
      (sigs(nmaj,nbins), &
       pe  (nmaj,nbins), &
       pin (nmaj,nbins))

    allocate                  &
      (sigex(nei,nmaj,nbins), &
       sigix(nei,nmaj,nbins))

    allocate                  &
      (siga(nei,nbins,nbins), &
       sec (nei,nbins,nbins))

    allocate(iimaxx(nbins))

    allocate           &
      (ww  (nei,nmaj), &
       ao  (nei,nmaj), &
       omeg(nei,nmaj), &
       anu (nei,nmaj), &
       bb  (nei,nmaj), &
       auto(nei,nmaj), &
       thi (nei,nmaj), &
       ak  (nei,nmaj), &
       aj  (nei,nmaj), &
       ts  (nei,nmaj), &
       ta  (nei,nmaj), &
       tb  (nei,nmaj), &
       gams(nei,nmaj), &
       gamb(nei,nmaj))

    allocate       &
     (alpha(jmax), &
      beta (jmax), &
      gamma(jmax), &
      psi  (jmax), &
      delz (jmax), &
      del2 (jmax), &
      dela (jmax), &
      delp (jmax), &
      delm (jmax), &
      dels (jmax), &
      den  (jmax))

  end subroutine cxglow_init
!--------------------------------------------------------------------------------

end module cxglow
