module readtgcm
  use netcdf
  implicit none
!
! Read tgcm (tiegcm or timegcm) history file, obtaining needed inputs 
! for glow.
!
  integer :: iyear,iday
  real :: ut
!
! Grid dimensions and coordinates:
!
  integer :: nlon,nlat,nlev,nilev          ! grid dimensions
  real,allocatable,dimension(:),save :: &  ! grid coordinates
    glon,glat,zlev,zilev
!
! 2d and 3d global fields at a single time (history):
!
  real,allocatable,dimension(:,:,:),save :: & ! 3d fields (nlon,nlat,nlev)
    tn,un,vn,o2,o1,no,ti,te,z,n2,zg,ne,n2d,n4s
  real,allocatable,dimension(:,:),save :: & ! 2d fields (nlon,nlat)
    cusp,drizzle,alfa,nflux,eflux
!
! Scalars (1d in time):
!
  real,allocatable,dimension(:),save :: &     ! scalars (ntime)
    f107d,f107a,e1,e2,h1,h2,alfac,ec,alfad,ed

  integer,parameter :: nf=19     ! number of fields to read
  character(len=8) :: fnames(nf) ! field names on file

  contains
!-----------------------------------------------------------------------
  subroutine read_tgcm(ncfile,mtime_request)
!
! Args:
    character(len=*),intent(in) :: ncfile
    integer,intent(in) :: mtime_request(3) ! requested model time
!
! Local:
    integer :: n,istat,ncid,k,i,ndims
    integer :: id,idunlim,idv_lon,idv_lat,idv_lev,idv_ilev,idv_time
    character(len=1024) :: msg
    character(len=NF90_MAX_NAME) :: varname
    integer :: ntime,itime,iyear1(1),iday1(1)
    real :: ut1(1),dlev
    real,allocatable,save :: f3d(:,:,:) ! global 3d field (lon,lat,lev)
    real,allocatable,save :: f2d(:,:)   ! global 2d field (lon,lat)
    real,allocatable,dimension(:,:,:),save :: o2mmr,o1mmr,n2mmr
    real,allocatable,save :: time(:)
    integer,allocatable,save :: mtime(:,:) ! (3,ntime)
    logical :: found(nf)
!
! These fields with these names must be on the file as 3d+time variables:
! Note that CUSP,DRIZZLE,ALFA,NFLUX,EFLUX are 2d (lon,lat).
!
    fnames = (/'TN      ','UN      ','VN      ','O2      ','O1      ',&
               'NO      ','TI      ','TE      ','Z       ','ZG      ',&
               'NE      ','N2      ','N2D     ','N4S     ','CUSP    ',&
               'DRIZZLE ','ALFA    ','NFLUX   ','EFLUX   '/)

    istat = nf90_open(ncfile,NF90_NOWRITE,ncid)
    if (istat /= NF90_NOERR) then
      write(msg,"('Error opening file ',a)") trim(ncfile)
      call handle_ncerr(istat,trim(msg),1)
    else
      write(6,"(/,'Opened file ',a)") trim(ncfile)
    endif
!
! Get number of times (histories) (length of unlimited variable):
!
    istat = nf90_inq_dimid(ncid,'time',idunlim)
    istat = nf90_inquire_dimension(ncid,idunlim,varname,ntime)
!   write(6,"('read_tgcm: ntime=',i4)") ntime
!
! Check for requested model time. Requested model time was initialized to
! -999,0,0. If it was not read by namelist (i.e., model day is still -999), 
! then use the first history on the file.
!
    allocate(mtime(3,ntime))
    istat = nf90_inq_varid(ncid,'mtime',id)
    istat = nf90_get_var(ncid,id,mtime,(/1,1/),(/3,ntime/)) 
    itime = 0
    if (mtime_request(1) == -999) then ! read first history on the file
      itime = 1
      write(6,"('Will read first history on the file (mtime ',3i4,').')") &
        mtime(:,1)
    else
      do i=1,ntime
        if (mtime(1,i)==mtime_request(1).and. &
            mtime(2,i)==mtime_request(2).and. &
            mtime(3,i)==mtime_request(3)) then
          itime = i
          write(6,"('Found requested mtime ',3i4,' (history ',i4,' on the file).')") &
            mtime_request,itime
          exit 
        endif
      enddo
    endif
    if (itime == 0) then
      write(6,"('>>> Could not find tgcm_mtime ',3i4,' on tgcm_ncfile ',a)") &
        mtime_request,trim(ncfile)
      write(6,"('There are ',i4,' histories on the file, with the following mtimes:')") ntime
      do i=1,ntime
        write(6,"('History ',i4,': mtime=',3i4)") i,mtime(:,i)
      enddo
      stop 'tgcm_mtime not found.'
    endif
!
! Get grid dimensions from the file:
!
    istat = nf90_inq_dimid(ncid,'lon',id)
    istat = nf90_inquire_dimension(ncid,id,varname,nlon)
    istat = nf90_inq_dimid(ncid,'lat',id)
    istat = nf90_inquire_dimension(ncid,id,varname,nlat)
    istat = nf90_inq_dimid(ncid,'lev',id)
    istat = nf90_inquire_dimension(ncid,id,varname,nlev)
    istat = nf90_inq_dimid(ncid,'ilev',id)
    istat = nf90_inquire_dimension(ncid,id,varname,nilev)
!
! Allocate coordinate variables:
!
    allocate(glon(nlon))
    allocate(glat(nlat))
    allocate(zlev(nlev))
    allocate(zilev(nilev))
    allocate(time(ntime))
!
! Read coordinate vars:
!
    istat = nf90_inq_varid(ncid,'lon',idv_lon)
    istat = nf90_get_var(ncid,idv_lon,glon,(/1/),(/nlon/))
    istat = nf90_inq_varid(ncid,'lat',idv_lat)
    istat = nf90_get_var(ncid,idv_lat,glat,(/1/),(/nlat/))
    istat = nf90_inq_varid(ncid,'lev',idv_lev)
    istat = nf90_get_var(ncid,idv_lev,zlev,(/1/),(/nlev/))
    istat = nf90_inq_varid(ncid,'ilev',idv_ilev)
    istat = nf90_get_var(ncid,idv_ilev,zilev,(/1/),(/nilev/))
    dlev = (zilev(nilev)-zilev(1))/(nilev-1)

    istat = nf90_inq_varid(ncid,'time',idv_time)
    istat = nf90_get_var(ncid,idv_time,time,(/1/),(/ntime/))

!   write(6,"('read_tgcm: nlon =',i4,' glon =',/,(8f10.3))") nlon,glon
!   write(6,"('read_tgcm: nlat =',i4,' glat =',/,(8f10.3))") nlat,glat
!   write(6,"('read_tgcm: nlev =',i4,' zlev =',/,(8f10.3))") nlev,zlev
!   write(6,"('read_tgcm: nilev=',i4,' zilev=',/,(8f10.3))") nilev,zilev
!   write(6,"('read_tgcm: nilev=',i4,' zilev(nilev)=',f9.3,' zilev(1)=',f9.3,' dlev=',f9.3)") &
!     nilev,zilev(nilev),zilev(1),dlev
    write(6,"('read_tgcm: ntime=',i4,' nlon=',i4,' nlat=',i4,' nlev=',i4)") ntime,nlon,nlat,nlev

    if (.not.allocated(f3d)) allocate(f3d(nlon,nlat,nlev)) ! 3d global read array
    if (.not.allocated(f2d)) allocate(f2d(nlon,nlat))      ! 2d global read array
!
! Allocate and read 1d time-dependent variables (or, this could be done
! on a per-history basis below, reading into scalars instead of arrays
! (like ut,iday,iyear, etc below).
!
!   f107d,f107a,e1,e2,h1,h2,alfac,ec,alfad,ed
!
    if (.not.allocated(f107d)) allocate(f107d(ntime))
    if (.not.allocated(f107a)) allocate(f107a(ntime))
    if (.not.allocated(e1))    allocate(e1(ntime))
    if (.not.allocated(e2))    allocate(e2(ntime))
    if (.not.allocated(h1))    allocate(h1(ntime))
    if (.not.allocated(h2))    allocate(h2(ntime))
    if (.not.allocated(alfac)) allocate(alfac(ntime))
    if (.not.allocated(alfad)) allocate(alfad(ntime))
    if (.not.allocated(ec))    allocate(ec(ntime))
    if (.not.allocated(ed))    allocate(ed(ntime))

    istat = nf90_inq_varid(ncid,'f107d',id)
    istat = nf90_get_var  (ncid,id,f107d)
    istat = nf90_inq_varid(ncid,'f107a',id)
    istat = nf90_get_var  (ncid,id,f107a)

    istat = nf90_inq_varid(ncid,'e1',id)
    istat = nf90_get_var  (ncid,id,e1)
    istat = nf90_inq_varid(ncid,'e2',id)
    istat = nf90_get_var  (ncid,id,e2)

    istat = nf90_inq_varid(ncid,'h1',id)
    istat = nf90_get_var  (ncid,id,h1)
    istat = nf90_inq_varid(ncid,'h2',id)
    istat = nf90_get_var  (ncid,id,h2)

    istat = nf90_inq_varid(ncid,'alfac',id)
    istat = nf90_get_var  (ncid,id,alfac)
    istat = nf90_inq_varid(ncid,'alfad',id)
    istat = nf90_get_var  (ncid,id,alfad)

    istat = nf90_inq_varid(ncid,'ec',id)
    istat = nf90_get_var  (ncid,id,ec)
    istat = nf90_inq_varid(ncid,'ed',id)
    istat = nf90_get_var  (ncid,id,ed)

!   write(6,"('read_tgcm: ntime=',i4,' f107d=',/,(8f9.2))") ntime,f107d
!   write(6,"('read_tgcm: ntime=',i4,' f107a=',/,(8f9.2))") ntime,f107a
!   write(6,"('read_tgcm: ntime=',i4,' e1=   ',/,(8f9.2))") ntime,e1
!   write(6,"('read_tgcm: ntime=',i4,' e2=   ',/,(8f9.2))") ntime,e2
!   write(6,"('read_tgcm: ntime=',i4,' h1=   ',/,(8f9.2))") ntime,h1
!   write(6,"('read_tgcm: ntime=',i4,' h2=   ',/,(8f9.2))") ntime,h2
!   write(6,"('read_tgcm: ntime=',i4,' alfac=',/,(8f9.2))") ntime,alfac
!   write(6,"('read_tgcm: ntime=',i4,' alfad=',/,(8f9.2))") ntime,alfad
!   write(6,"('read_tgcm: ntime=',i4,' ec=   ',/,(8f9.2))") ntime,ec
!   write(6,"('read_tgcm: ntime=',i4,' ed=   ',/,(8f9.2))") ntime,ed
!
! Read itime history (model time mtime(:,itime)):
!
      istat = nf90_inq_varid(ncid,'year',id)
      istat = nf90_get_var(ncid,id,iyear1,(/itime/),(/1/))
      iyear = iyear1(1)

      istat = nf90_inq_varid(ncid,'day',id)
      istat = nf90_get_var(ncid,id,iday1,(/itime/),(/1/))
      iday = iday1(1)

      istat = nf90_inq_varid(ncid,'ut',id)
      istat = nf90_get_var(ncid,id,ut1,(/itime/),(/1/))
      ut = ut1(1)

      write(6,"(/,'read_tgcm: itime=',i4,' ntime=',i4,' iyear=',i5,' iday=',i4,' ut=',f8.2,' mtime=',3i4)") &
        itime,ntime,iyear,iday,ut,mtime(:,itime)
!
! Read vars at current history/time:
!
      found(:) = .false.
      do n=1,nf
        istat = nf90_inq_varid(ncid,fnames(n),id)
        if (istat /= NF90_NOERR) then
          write(6,"('>>> WARNING: did not find field ',a,' on the tgcm file')") &
            trim(fnames(n))
          cycle
        endif
        istat = nf90_inquire_variable(ncid,id,ndims=ndims)
!
! 3d var+time: assume (nlon,nlat,nlev):
!
        if (ndims==4) then
          f3d = 0.
          istat = nf90_get_var(ncid,id,f3d,(/1,1,1,itime/),(/nlon,nlat,nlev,1/))
          if (istat /= NF90_NOERR) then
            write(6,"('>>> Error reading 3d var ')") fnames(n)
            call handle_ncerr(istat,'Error from nf90_get_var',1)
          endif
!         write(6,"('Read 3d var: itime=',i4,' n=',i4,' field ',a,' id=',i4,' min,max=',2es12.4)") &
!           itime,n,trim(fnames(n)),id,minval(f3d(:,:,1:nlev-1)),maxval(f3d(:,:,1:nlev-1))
!
! 2d var+time: assume (nlon,nlat):
!
        elseif (ndims==3) then
          f2d = 0.
          istat = nf90_get_var(ncid,id,f2d,(/1,1,itime/),(/nlon,nlat,1/))
          if (istat /= NF90_NOERR) then
            write(6,"('>>> Error reading 2d var ')") fnames(n)
            call handle_ncerr(istat,'Error from nf90_get_var',1)
          endif
!         write(6,"('Read 2d var: time=',i4,' n=',i4,' field ',a,' id=',i4,' min,max=',2es12.4)") &
!           itime,n,trim(fnames(n)),id,minval(f2d(:,:)),maxval(f2d(:,:))
        endif ! 2d or 3d var
!
! Transfer field to module data:
! Top level may be missing in tn,un,vn,ti,te, so set these to values at nlev-1 for now:
!
        select case(trim(fnames(n)))
          case('TN')
            if (.not.allocated(tn)) allocate(tn(nlon,nlat,nlev))
            tn = f3d
            tn(:,:,nlev) = tn(:,:,nlev-1) 
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(tn),maxval(tn)
            found(n) = .true. 
          case('UN')
            if (.not.allocated(un)) allocate(un(nlon,nlat,nlev))
            un = f3d
            un(:,:,nlev) = un(:,:,nlev-1) 
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(un),maxval(un)
            found(n) = .true. 
          case('VN')
            if (.not.allocated(vn)) allocate(vn(nlon,nlat,nlev))
            vn = f3d
            vn(:,:,nlev) = vn(:,:,nlev-1) 
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(vn),maxval(vn)
            found(n) = .true. 
          case('O2')
            if (.not.allocated(o2)) allocate(o2(nlon,nlat,nlev))
            o2 = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(o2),maxval(o2)
            found(n) = .true. 
          case('O1')
            if (.not.allocated(o1)) allocate(o1(nlon,nlat,nlev))
            o1 = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(o1),maxval(o1)
            found(n) = .true. 
          case('N2')
            if (.not.allocated(n2)) allocate(n2(nlon,nlat,nlev))
            n2 = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(n2),maxval(n2)
            found(n) = .true. 
          case('NO')
            if (.not.allocated(no)) allocate(no(nlon,nlat,nlev))
            no = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(no),maxval(no)
            found(n) = .true. 
          case('TI')
            if (.not.allocated(ti)) allocate(ti(nlon,nlat,nlev))
            ti = f3d
            ti(:,:,nlev) = ti(:,:,nlev-1) 
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(ti),maxval(ti)
            found(n) = .true. 
          case('TE')
            if (.not.allocated(te)) allocate(te(nlon,nlat,nlev))
            te = f3d
            te(:,:,nlev) = te(:,:,nlev-1) 
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(te),maxval(te)
            found(n) = .true. 
          case('NE')
            if (.not.allocated(ne)) allocate(ne(nlon,nlat,nlev))
            ne = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(ne),maxval(ne)
            found(n) = .true. 
          case('Z')
            if (.not.allocated(z)) allocate(z(nlon,nlat,nlev))
            z = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(z),maxval(z)
            found(n) = .true. 
          case('ZG')
            if (.not.allocated(zg)) allocate(zg(nlon,nlat,nlev))
            zg = f3d*1.e-5 ! cm->km
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(zg),maxval(zg)
            found(n) = .true. 
          case('N2D')
            if (.not.allocated(n2d)) allocate(n2d(nlon,nlat,nlev))
            n2d = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(n2d),maxval(n2d)
            found(n) = .true. 
          case('N4S')
            if (.not.allocated(n4s)) allocate(n4s(nlon,nlat,nlev))
            n4s = f3d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(n4s),maxval(n4s)
            found(n) = .true. 
!
! 2d fields (nlon,nlat):
!
          case('CUSP')
            if (.not.allocated(cusp)) allocate(cusp(nlon,nlat))
            cusp = f2d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(cusp),maxval(cusp)
            found(n) = .true. 
          case('DRIZZLE')
            if (.not.allocated(drizzle)) allocate(drizzle(nlon,nlat))
            drizzle = f2d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(drizzle),maxval(drizzle)
            found(n) = .true. 
          case('ALFA')
            if (.not.allocated(alfa)) allocate(alfa(nlon,nlat))
            alfa = f2d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(alfa),maxval(alfa)
            found(n) = .true. 
          case('NFLUX')
            if (.not.allocated(nflux)) allocate(nflux(nlon,nlat))
            nflux = f2d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(nflux),maxval(nflux)
            found(n) = .true. 
          case('EFLUX')
            if (.not.allocated(eflux)) allocate(eflux(nlon,nlat))
            eflux = f2d
            write(6,"('read_tgcm: Field ',a,' min,max=',2es12.4)") &
              fnames(n),minval(eflux),maxval(eflux)
            found(n) = .true. 
!
! Unknown field (this should not happen):
!
          case default
            write(6,"('>>> read_tgcm: unknown field: ',a)") trim(fnames(n))
            stop 'read tgcm 3d fields'
        end select ! fnames(n)
      enddo ! n=1,nf
!
! Define n2 mmr from o2,o if necessary:
!
      if (.not.found(findx('N2'))) then
        if (.not.allocated(N2)) allocate(n2(nlon,nlat,nlev))
        if (found(findx('O1')).and.found(findx('O2'))) then
          n2 = 1.-o2-o1
          write(6,"('read_tgcm: Field N2 (1-O2-O) min,max=',2es12.4)") &
            minval(n2),maxval(n2)
        else
          write(6,"('>>> WARNING read_tgcm: O2 or O apparently not read: am setting n2=0.')")
          n2 = 0.
        endif
        found(findx('N2'))=.true.
      endif ! N2 not found
!
! If ZG not on the file, calculate it using sub calczg (from tgcm addiag.F):
!
      if (.not.found(findx('ZG'))) then
        if (.not.found(findx('TN')).or..not.found(findx('O2')).or.&
            .not.found(findx('O1')).or..not.found(findx('N2')).or.&
            .not.found(findx('Z' ))) then
          write(6,"('>>> FATAL: ZG not found on the file, and at least one of')")
          write(6,"('    tn,o2,o,n2,z also not found on the file.')")
          stop 'read_tgcm'
        endif
        if (.not.allocated(zg)) allocate(zg(nlon,nlat,nlev))
        call calczg(tn,o2,o1,n2,z,zg,glat,nlon,nlat,nlev,dlev)
        zg = zg*1.e-5 ! cm to km
        write(6,"('Field ZG (from calczg) min,max (km)=',2es12.4)") &
          minval(zg),maxval(zg)
        found(findx('ZG'))=.true.
      endif ! zg not on the file
!
! If N2D was not found, set it to zero:
!
      if (.not.found(findx('N2D'))) then
        if (.not.allocated(n2d)) allocate(n2d(nlon,nlat,nlev))
        n2d = 0.
        write(6,"('read_tgcm: N2D not found on tgcm history, so set it to zero.')")
        found(findx('N2D'))=.true.
      endif
!
! Auroral fields not found can be set to zero:
!
      if (.not.found(findx('CUSP'))) then
        if (.not.allocated(cusp)) allocate(cusp(nlon,nlat))
        cusp = 0.
        write(6,"('read_tgcm: CUSP not found on tgcm history, so set it to zero.')")
        found(findx('CUSP'))=.true.
      endif
      if (.not.found(findx('DRIZZLE'))) then
        if (.not.allocated(drizzle)) allocate(drizzle(nlon,nlat))
        drizzle = 0.
        write(6,"('read_tgcm: DRIZZLE not found on tgcm history, so set it to zero.')")
        found(findx('DRIZZLE'))=.true.
      endif
      if (.not.found(findx('ALFA'))) then
        if (.not.allocated(alfa)) allocate(alfa(nlon,nlat))
        alfa = 0.
        write(6,"('read_tgcm: ALFA not found on tgcm history, so set it to zero.')")
        found(findx('ALFA'))=.true.
      endif
      if (.not.found(findx('NFLUX'))) then
        if (.not.allocated(nflux)) allocate(nflux(nlon,nlat))
        nflux = 0.
        write(6,"('read_tgcm: NFLUX not found on tgcm history, so set it to zero.')")
        found(findx('NFLUX'))=.true.
      endif
      if (.not.found(findx('EFLUX'))) then
        if (.not.allocated(eflux)) allocate(eflux(nlon,nlat))
        eflux = 0.
        write(6,"('read_tgcm: EFLUX not found on tgcm history, so set it to zero.')")
        found(findx('EFLUX'))=.true.
      endif
!
! If any other fields were not found on the history, stop w/ fatal error:
!
      do n=1,nf
        if (.not.found(n)) then
          write(6,"(/,'>>> FATAL: field ',a,' was not found on tgcm file ',a)") &
            fnames(n),trim(ncfile)
          stop 'read_tgcm'
        endif
      enddo
!
! Convert species from mmr to cm3 (still in the time/history loop):
!
      if (.not.allocated(o2mmr)) allocate(o2mmr(nlon,nlat,nlev))
      if (.not.allocated(o1mmr)) allocate(o1mmr(nlon,nlat,nlev))
      if (.not.allocated(n2mmr)) allocate(n2mmr(nlon,nlat,nlev))
      o2mmr = o2
      o1mmr = o1
      n2mmr = n2

      call denconv(o2, 32.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)
      call denconv(o1, 16.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)
      call denconv(n2, 28.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)
      call denconv(no, 28.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)
      call denconv(n2d,14.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)
      call denconv(n4s,14.,tn,o2mmr,o1mmr,n2mmr,zlev,nlon,nlat,nlev)

      write(6,"('After denconv: O2 (cm3) min,max=',2es12.4)") minval(o2),maxval(o2)
      write(6,"('After denconv: O1 (cm3) min,max=',2es12.4)") minval(o1),maxval(o1)
      write(6,"('After denconv: N2 (cm3) min,max=',2es12.4)") minval(n2),maxval(n2)
      write(6,"('After denconv: NO (cm3) min,max=',2es12.4)") minval(no),maxval(no)
      write(6,"('After denconv: N2D(cm3) min,max=',2es12.4)") minval(n2d),maxval(n2d)
      write(6,"('After denconv: N4S(cm3) min,max=',2es12.4)") minval(n4s),maxval(n4s)
!
! Close the file:
!
    istat = nf90_close(ncid)
  end subroutine read_tgcm
!-----------------------------------------------------------------------
      subroutine calczg(tn,o2,o1,n2,z,zg,glat,nlon,nlat,nlev,dlev)

!     use params_module,only: dz,glat
!     use init_module,only: istep
!     use cons_module,only: boltz,avo
!     use addfld_module,only: addfld
!
! Given geopotential z (calculated with the model constant gravity),
!   calculate geopotential zg with varying gravity. This is taken from
!   tgcmproc_f90, routines calchts and glatf in proclat.F.
!
! Args:
      integer,intent(in) :: nlon,nlat,nlev
      real,dimension(nlon,nlat,nlev),intent(in):: &
        tn, & ! neutral temperature (deg K)
        o2, & ! molecular oxygen (mmr)
        o1, & ! atomic oxygen (mmr)
        n2, & ! molecular nitrogen (mmr)
        z     ! geopotential calculated with constant gravity (from addiag)
      real,intent(in) :: glat(nlat),dlev
      real,dimension(nlon,nlat,nlev),intent(out) :: &
        zg  ! output geopotential calculated with varying gravity
!
! Local:
      integer :: i,j,k
      real :: g0,r0,c2
      real,dimension(nlev) :: xmas,g
      real,parameter :: dgtr = 1.74533E-2
      real,parameter :: avo = 6.023e23
      real,parameter :: boltz = 1.38E-16
!
! Latitude scan:
      do j=1,nlat
        c2 = cos(2.*dgtr*glat(j))
        g0 = 980.616*(1.-.0026373*c2)
        r0 = 2.*g0/(3.085462e-6 + 2.27e-9*c2) ! effective earth radius
!
! Longitude scan:
        do i=1,nlon
          g(1)=g0*(r0/(r0+0.5*(z(i,j,1)+z(i,j,2))))**2
          xmas(:) = 1./(o1(i,j,:)/16.+o2(i,j,:)/32.+ &
                        n2(i,j,:)/28.)/avo
!
! Levels:
          zg(i,j,1) = z(i,j,1)
          do k=2,nlev-1
            zg(i,j,k) = zg(i,j,k-1) + boltz*dlev*tn(i,j,k-1) / &
              (xmas(k-1)*g(k-1))
            g(k)=g0*(r0/(r0+0.5*(zg(i,j,k)+z(i,j,k+1))))**2
          enddo ! k=lev0+1,lev1-1
          zg(i,j,nlev) = 1.5*zg(i,j,nlev-1)-0.5*zg(i,j,nlev-2)
        enddo ! i=lon0,lon1
      enddo ! j=1,nlat
      end subroutine calczg
!-----------------------------------------------------------------------
  subroutine denconv(f,wt,tn,o2,o1,n2,zlev,nlon,nlat,nlev)
!
! Convert density input field f from mmr to cm3.
! Inputs o2,o1,n2 are in mmr.
! zlev is the vertical coordinate (zp).
!
  implicit none
!
! Args:
  integer,intent(in) :: nlon,nlat,nlev
  real,intent(inout) :: f(nlon,nlat,nlev)
  real,dimension(nlon,nlat,nlev),intent(in) :: tn,o2,o1,n2 ! mmr
  real,intent(in) :: wt
!
! Assume levels coordinate is increasing from bottom to top, regular dz)
  real,intent(in) :: zlev(nlev) ! levels coord (tgcm zp)
!
! Local:
  integer :: i,j,k
  real :: zp,dlev,mbar,pkt
  real,parameter :: p0=5.e-4 ! microbars
  real,parameter :: boltz=1.3805e-16

  dlev = (zlev(nlev)-zlev(1))/(nlev-1)
  do k=1,nlev
    zp = zlev(1)+(k-1)*dlev
    do j=1,nlat
      do i=1,nlon
        mbar = 1./(o2(i,j,k)/32.+o1(i,j,k)/16.+n2(i,j,k)/28.)
        pkt = p0*exp(-zp)/(boltz*tn(i,j,k))
        f(i,j,k) = f(i,j,k) * pkt * mbar / wt
      enddo
    enddo
  enddo 
  end subroutine denconv
!-----------------------------------------------------------------------
  integer function findx(name)
!
! Find index in fnames(nf) that matches input name.
!
  character(len=*),intent(in) :: name
  integer :: n    
  findx = 0
  do n=1,nf
    if (trim(fnames(n))==name) then
      findx = n
      exit
    endif
  enddo
  if (findx==0) then
    write(6,"('>>> findx: could not find index to field ',a)") name
    stop 'findx'
  endif
  end function findx
!-----------------------------------------------------------------------
  subroutine handle_ncerr(istat,msg,ifatal)
!
! Handle a netcdf lib error:
!
    integer,intent(in) :: istat,ifatal
    character(len=*),intent(in) :: msg
!
    write(6,"(/72('-'))")
    write(6,"('>>> Error from netcdf library:')")
    write(6,"(a)") trim(msg)
    write(6,"('istat=',i5)") istat
    write(6,"(a)") nf90_strerror(istat)
    write(6,"(72('-')/)")
    if (ifatal > 0) stop('Fatal netcdf error')
  end subroutine handle_ncerr
!-----------------------------------------------------------------------
end module readtgcm
