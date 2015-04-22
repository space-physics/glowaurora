program glow_drv

! Main driver for the GLOW model.
! Uses TIE-GCM history files or MSIS/IRI for input.

! Stan Solomon and Ben Foster, 1/15

! For definitions of shared variables, see subroutine GLOW and module CGLOW.

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

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
! NF      number of types of auroral fluxes (obsolete)

  use readtgcm,only: read_tgcm    ! subroutine to read tgcm history files
  use readtgcm,only: nlon_tgcm=>nlon, glon_tgcm=>glon
  use readtgcm,only: nlat_tgcm=>nlat, glat_tgcm=>glat
  use readtgcm,only: nlev_tgcm=>nlev
  use readtgcm,only: iyear_tgcm=>iyear,iday_tgcm=>iday,ut_tgcm=>ut
  use readtgcm,only: f107_tgcm=>f107d,f107a_tgcm=>f107a
  use readtgcm,only: alfacusp=>alfac,ecusp=>ec,alfadrizzle=>alfad,edrizzle=>ed
  use readtgcm,only: eflux,nflux,alfa,drizzle,cusp
  use readtgcm,only: zg,tn,un,vn,o2,o1,n2,n4s,n2d,no,ti,te,ne

  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst,nf
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac,hlybr,fexvir,hlya,heiew
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,ecalc,zxden,zeta,zceta
  use cglow,only: cglow_init
  use cxglow,only: cxglow_init

  use output,only: output_init    ! subroutine to allocate output arrays
  use output,only: create_ncfile  ! subroutine to create netcdf output file
  use output,only: write_ncfile   ! subroutine to write data to netcdf output file
  use output,only: nlon,nlat,nlev ! grid dimensions (set here in glow_drv)
  use output,only: lon,lat,lev    ! grid coordinates (set here in glow_drv)
  use output,only: zzz,ao,ao2,an2,ano,an4s,an2d,atn,ati,ate,ane,aun,avn,nec,xden,eta ! output arrays

  implicit none

  character(len=1024) :: &
    tgcm_ncfile,         &    ! path to tgcm history file (tiegcm or timegcm)
    glow_ncfile,         &    ! path to netcf glow output file
    iri90_dir                 ! directory containing iri data files

  real,allocatable :: z(:)            ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:),zvn(:)   ! winds on glow grid
  real,allocatable :: outf(:,:)       ! iri output (11,jmax)
  real :: rz12,stl,glati,glongi,fmono,emono
  real :: d(8), t(2), sw(25), oarr(30)
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,itail
  logical :: tgcm, loop_lons, first, jf(12)
  integer :: tgcm_mtime(3) ! tgcm model time day,hour,minute
  data sw/25*1./, loop_lons/.true./, first/.true./

  namelist /glow_input/ &
    idate,ut,glati,glongi,f107a,f107,f107p,ap, &
    iscale,jlocal,kchem,xuvfac,ef,ec,itail,fmono,emono, &
    tgcm_ncfile,iri90_dir,loop_lons,jmax,glow_ncfile, &
    tgcm_mtime

! Execute:
!
! Set obsolete parameters to zero:
!
  hlybr = 0.   ! obsolete
  fexvir = 0.  ! obsolete
  hlya = 0.    ! obsolete
  heiew = 0.   ! obsolete
!
! Initialize tgcm_ncfile; if not provided, will do an MSIS/IRI run.
! If tgcm_mtime is not read from namelist, model day of -999 
! will flag read_tgcm to read the first history on the file.
!
  tgcm_ncfile = ' '
  tgcm_mtime = (/-999,0,0/) 
!
! Namelist read from input file:
!
  read (5,nml=glow_input)
!
! Read tgcm history file (tiegcm or timegcm), if provided:
!
  if (len_trim(tgcm_ncfile) > 0) then
    tgcm = .true.
    write(6,"('glow_drv: tgcm_ncfile = ',a)") trim(tgcm_ncfile)
    write(6,"('glow_drv: tgcm_mtime  = ',3i4)") tgcm_mtime
    call read_tgcm(tgcm_ncfile,tgcm_mtime)
    jmax=68
    nlev = jmax
    nlat = nlat_tgcm
    nlon = nlon_tgcm
  else
    tgcm = .false.
    jmax = 86
    nlev = jmax
    nlat = 36
    nlon = 72
    write(6,"('glow_drv: tgcm_ncfile not provided, will use MSIS/IRI')")
  endif

  if (loop_lons == .false.) then
    nlon = 1
    write(6,"('glow_drv: loop_lons = .false. --> will not loop over longitudes.')")
  endif
!
! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(outf(11,jmax))
!
! Allocate arrays in other modules (formerly in common blocks):
!
  call cglow_init
  call cxglow_init
  call output_init
!
! Set glow grid coordinates:
!
    do k=1,nlev
      lev(k)=float(k)/4.-10.25
    enddo
  if (tgcm) then
    lon(:) = glon_tgcm(:)
    lat(:) = glat_tgcm(:)
  else
    if (loop_lons == .true.) then
      do i=1,nlon
        lon(i) = float(i)*5. - 185. ! note this is tgcm 5-deg grid 
      enddo
    else
      lon(1) = glongi ! from namelist
    endif
    do j=1,nlat
      lat(j)=j*5.-92.5
    enddo
  endif
!
! Set up electron energy grid:
!
  call egrid (ener, del, nbins)
!
! Set input parameters if this is a tgcm run (otherwise use namelist inputs):
!
  if (tgcm) then
    idate=iyear_tgcm*1000+iday_tgcm
    ut   =ut_tgcm*3600.
    f107 =f107_tgcm(1)
    f107p=f107_tgcm(1)
    f107a=f107a_tgcm(1)
  endif

  write(6,"('glow_drv: idate=',i7,' ut=',f7.1)") idate,ut
  write(6,"('glow_drv: f107a=',f5.1,' f107=',f5.1,' f107p=',f5.1,' ap=',f5.1)") &
    f107a,f107,f107p,ap
!
! Loop over longitude:
!
  do i=1,nlon
    glong = lon(i)
!
! Loop over latitude:
!
    do l=1,nlat
      glat = lat(l)
!
! Calculate local solar time:
!
      stl = ut/3600. + glong/15.
      if (stl < 0.) stl = stl + 24.
      if (stl >= 24.) stl = stl - 24.
!
! If this is a tgcm run, use altitude grid and inputs from tgcm history:
!
      if (tgcm) then
        call tzgrid(i,l,jmax,z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte)
!
! Put auroral electron flux from tgcm history into phitop array:
!
        phitop(:) = 0.
        ef=eflux(i,l)
        ec=alfa(i,l)*1000.    ! keV to eV
        if (ef>.001 .and. ec>1.) call maxt(ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
!
! If this is not a tgcm run, use default grid and MSIS/NOEM/IRI inputs:
!
      else

        call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                   z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Put auroral electron flux from namelist input into phitop array:
!
        phitop(:) = 0.
        if (ef>.001 .and. ec>1.) call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
      endif 
!
! Fill altitude array, converting back to cm:
! (note that zg from tgcm was converted to km in read_tgcm.F90)
!
      zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission
! rates, and vertical column brightnesses:
!
      call glow
!
! Collect arrays for output to netCDF:
!
      do j=1,jmax
        zzz(i,l,j)  = zz(j)
        ao(i,l,j)   = zo(j)
        ao2(i,l,j)  = zo2(j)
        an2(i,l,j)  = zn2(j)
        ano(i,l,j)  = zno(j)
        an4s(i,l,j) = zns(j)
        an2d(i,l,j) = znd(j)
        atn(i,l,j)  = ztn(j)
        ati(i,l,j)  = zti(j)
        ate(i,l,j)  = zte(j)
        aun(i,l,j)  = zun(j)
        avn(i,l,j)  = zvn(j)
        ane(i,l,j)  = ze(j)
        nec(i,l,j)  = ecalc(j)
        do ix=1,nex
          xden(i,l,j,ix) = zxden(ix,j)
        enddo
        do ix=1,nw
          eta(i,l,j,ix) = zeta(ix,j)
        enddo
      enddo

    enddo    ! latitude loop

  enddo   ! longitude loop

!
! Output section:
!
! Create and define netcdf output file, if this is the first history processed:
!
  if (first == .true.) call create_ncfile(glow_ncfile,tgcm_ncfile)
  first = .false.
!
! Write global arrays to netCDF file:
!
  call write_ncfile(glow_ncfile)

end program glow_drv

!----------------------------------------------------------------------

! Subroutine TZGRID maps fields from TIE-GCM or TIME-GCM history onto
! GLOW altitude grid, interpolating and/or extrapolating as necessary.

! GLOW altitude grid extends from lev=ln(P0/P) -10 to +6.75, in intervals of H/4.

! Lower boundary extrapolation for TIE-GCM from lev =-10 to -7 assumes H = 6 km

! low-res TIE-GCM: 29 levels; extrapolate lower boundary, use both midpoints and interfaces
! high-res TIE-GCM: 57 levels; extrapolate lower boundary, use interfaces
! low-res TIME-GCM: 49 levels; use both midpoints and interfaces
! high-res TIME-GCM: 97 levels; use interfaces

subroutine tzgrid(i,l,jmax,z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte)

  use readtgcm,only: nlev,zg,tn,un,vn,o2,o1,n2,n4s,n2d,no,ti,te,ne

  implicit none

  integer,intent(in) :: i,l,jmax
  real,intent(out) :: z(jmax),zo(jmax),zo2(jmax),zn2(jmax),zns(jmax),znd(jmax), &
                      zno(jmax),ztn(jmax),zti(jmax),zte(jmax),zun(jmax),zvn(jmax),ze(jmax)
  integer :: j,jj

  if (nlev/=29 .and. nlev/=57 .and. nlev/=49 .and. nlev/=97) then
    write(6,"('zgrid: unknown NLEV = ',i5)") nlev
    stop 'zgrid'
  endif

  if (nlev == 29) then                              ! low-res TIE-GCM
    do j=1,jmax
      jj=(j-11)/2
      if (j <= 13) then                             ! at or below lower boundary
        z(j)   = zg(i,l,1)-float(13-j)*6./4.        ! extrapolate
        zo(j)  = o1(i,l,1)*exp(-float(14-j)/4.)
        zo2(j) = o2(i,l,1)*exp(float(14-j)/4.)
        zn2(j) = n2(i,l,1)*exp(float(14-j)/4.)
        zns(j) = n4s(i,l,1)
        znd(j) = n2d(i,l,1)
        zno(j) = no(i,l,1)
        ztn(j) = tn(i,l,1)
        zti(j) = ti(i,l,1)
        zte(j) = te(i,l,1)
        zun(j) = un(i,l,1)
        zvn(j) = vn(i,l,1)
        ze(j)  = ne(i,l,1)
      else
        if (j/2*2 == j) then                        ! at midpoint
          z(j)   = (zg(i,l,jj)+zg(i,l,jj+1))/2      ! interpolate zg to midpoint
          zo(j)  = o1(i,l,jj)                       ! use fields at midpoint
          zo2(j) = o2(i,l,jj)
          zn2(j) = n2(i,l,jj)
          zns(j) = n4s(i,l,jj)
          znd(j) = n2d(i,l,jj)
          zno(j) = no(i,l,jj)
          ztn(j) = tn(i,l,jj)
          zti(j) = ti(i,l,jj)
          zte(j) = te(i,l,jj)
          zun(j) = un(i,l,jj)
          zvn(j) = vn(i,l,jj)
          ze(j)  = sqrt(ne(i,l,jj)*ne(i,l,jj+1))
        else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
        endif
      endif
    enddo 
    z(jmax) = z(jmax-1) + (z(jmax-1)-z(jmax-2))     !patch top level
    ze(jmax) = ze(jmax-1)**2 / ze(jmax-2)
  endif

  if (nlev == 57) then                              ! high-res TIE-GCM
    do j=1,jmax
      jj=j-12
      if (j <= 13) then                             ! at or below lower boundary
        z(j)   = zg(i,l,1)-float(13-j)*6./4.        ! extrapolate
        zo(j)  = o1(i,l,1)*exp( 0.125-float(14-j)/4.)
        zo2(j) = o2(i,l,1)*exp(-0.125+float(14-j)/4.)
        zn2(j) = n2(i,l,1)*exp(-0.125+float(14-j)/4.)
        zns(j) = n4s(i,l,1)
        znd(j) = n2d(i,l,1)
        zno(j) = no(i,l,1)
        ztn(j) = tn(i,l,1)
        zti(j) = ti(i,l,1)
        zte(j) = te(i,l,1)
        zun(j) = un(i,l,1)
        zvn(j) = vn(i,l,1)
        ze(j)  = ne(i,l,1)
      else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
      endif
    enddo 
  endif

  if (nlev == 49) then                              ! low-res TIME-GCM
    do j=1,jmax
        jj=(j+29)/2
        if (j/2*2 == j) then                        ! at midpoint
          z(j)   = (zg(i,l,jj)+zg(i,l,jj+1))/2      ! interpolate zg to midpoint
          zo(j)  = o1(i,l,jj)                       ! use fields at midpoint
          zo2(j) = o2(i,l,jj)
          zn2(j) = n2(i,l,jj)
          zns(j) = n4s(i,l,jj)
          znd(j) = n2d(i,l,jj)
          zno(j) = no(i,l,jj)
          ztn(j) = tn(i,l,jj)
          zti(j) = ti(i,l,jj)
          zte(j) = te(i,l,jj)
          zun(j) = un(i,l,jj)
          zvn(j) = vn(i,l,jj)
          ze(j)  = sqrt(ne(i,l,jj)*ne(i,l,jj+1))
        else
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
        endif
    enddo 
    z(jmax) = z(jmax-1) + (z(jmax-1)-z(jmax-2))     !patch top level
    ze(jmax) = ze(jmax-1)**2 / + ze(jmax-2)
  endif

  if (nlev == 97) then                              ! high-res TIME-GCM
    do j=1,jmax
          jj=j+28
          z(j)   = zg(i,l,jj)                       ! at interface
          zo(j)  = sqrt(o1(i,l,jj-1)*o1(i,l,jj))    ! interpolate fields to interface
          zo2(j) = sqrt(o2(i,l,jj-1)*o2(i,l,jj))
          zn2(j) = sqrt(n2(i,l,jj-1)*n2(i,l,jj))
          zns(j) = sqrt(n4s(i,l,jj-1)*n4s(i,l,jj))
          znd(j) = sqrt(n2d(i,l,jj-1)*n2d(i,l,jj))
          zno(j) = sqrt(no(i,l,jj-1)*no(i,l,jj))
          ztn(j) = (tn(i,l,jj-1)+tn(i,l,jj))/2.
          zti(j) = (ti(i,l,jj-1)+ti(i,l,jj))/2.
          zte(j) = (te(i,l,jj-1)+te(i,l,jj))/2.
          zun(j) = (un(i,l,jj-1)+un(i,l,jj))/2.
          zvn(j) = (vn(i,l,jj-1)+vn(i,l,jj))/2.
          ze(j)  = ne(i,l,jj)
    enddo 
  endif

end subroutine tzgrid

!----------------------------------------------------------------------

! Subroutine MZGRID gets fields from eMpirical models on default GLOW altitude grid.

! Neutral densities from NRLMSISE-00 model (a.k.a. MSIS2K).
! Nitric oxide densities from NOEM model (via snoem.f and snoemint.f).
! Electron densities from IRI-90.

subroutine mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                   z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)

  implicit none

  integer,intent(in) :: jmax,nex,idate
  real,intent(in) :: ut,glat,glong,stl,f107a,f107,f107p,ap
  character(len=1024),intent(in) :: iri90_dir
  real,intent(out) :: z(jmax),zo(jmax),zo2(jmax),zn2(jmax),zns(jmax),znd(jmax), &
       zno(jmax),ztn(jmax),zti(jmax),zte(jmax),zun(jmax),zvn(jmax),ze(jmax),zxden(nex,jmax)

  integer :: j,ijf,jmag,iday,mmdd
  real :: rz12, d(8), t(2), sw(25), oarr(30)
  logical :: jf(12)
  real,allocatable :: outf(:,:)              ! iri output (11,jmax)
  data sw/25*1./

  if (jmax /= 86) then
    write(6,"('mzgrid: unknown JMAX = ',i5)") jmax
    stop 'mzgrid'
  endif

  allocate(outf(11,jmax))

!
! Set default altitudes:
!
      z = (/ 80.0,81.5,83.0,84.5,86.0,87.5,89.0,90.5,92.0,93.5, &
             95.0,97.5,99.0,100.5,102.0,103.5,105.0,107.5,109.,111., &
             113.,115.,117.5,120.,122.5,125.,128.,131.,135.,139., &
             143.,148.,153.,168.,174.,180.,186.,193.,200.,207., &
             215.,223.,231.,239.,247.,255.,264.,273.,282.,291., &
             300.,310.,320.,330.,340.,350.,360.,370.,380.,390., &
             400.,410.,420.,430.,440.,450.,460.,470.,480.,490., &
             500.,510.,520.,530.,540.,550.,560.,570.,580.,590., &
             600.,610.,620.,630.,640.,650. /)
!
! Call MSIS-2K to get neutral densities and temperature:
!
        call tselec(sw)

        do j=1,jmax ! levels
          call gtd7(idate,ut,z(j),glat,glong,stl,f107a,f107p,ap,48,d,t)
          zo(j) = d(2)
          zn2(j) = d(3)
          zo2(j) = d(4)
          zns(j) = d(8)
          ztn(j) = t(2)
          znd(j)  = 0.
        enddo
!
! Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical Model (NOEM):
!
        call snoemint(idate,glat,glong,f107,ap,jmax,z,ztn,zno)
!
! Call International Reference Ionosphere-1990 subroutine to get
! electron density, electron temperature, and ion temperature:
! The directory iri90_dir is the location of the ccirnn.asc and ursinn.asc files.
!
        jf(:) = .true.
!       jf(12) = .false.
        jf(5) = .false.
        jmag = 0
        rz12 = -f107a
        iday = idate - idate/1000*1000
        mmdd = -iday
        outf = 0.

        call iri90(jf,jmag,glat,glong,rz12,mmdd,stl,z,jmax,trim(iri90_dir),outf,oarr)

        do j=1,jmax
          ze(j) = outf(1,j) / 1.E6
          if (ze(j) .lt. 100.) ze(j) = 100.
          zti(j) = outf(3,j)
          if (zti(j) .lt. ztn(j)) zti(j) = ztn(j)
          zte(j) = outf(4,j)
          if (zte(j) .lt. ztn(j)) zte(j) = ztn(j)
          zxden(3,j) = ze(j) * outf(5,j)/100.
          zxden(6,j) = ze(j) * outf(8,j)/100.
          zxden(7,j) = ze(j) * outf(9,j)/100.
        enddo
!
! Until implementation of an empirical neutral wind model, winds are set to zero:
!
        zun(:) = 0.
        zvn(:) = 0.

end subroutine mzgrid
