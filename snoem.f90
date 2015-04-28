! Subroutine SNOEM calculates nitric oxide zonal mean altitude profile
! as function of magnetic latitude for specified day of year, Kp, and F10.7.

! The NOEM empirical model is based on data from the SNOE ultraviolet
! spectrometer during 1998-2000, using empirical orthogonal function analysis.
! Altitude range is from 100 to 150 km.

! Marsh et al., JGR, 109, A07301, doi:10.1029/2003JA010199, 2004.

! Adapted by Stan Solomon, 5/14, from IDL and F90 code supplied by Dan Marsh. 
module snoemmod
    use machprec
    implicit none
    private
    public :: snoem
contains
  subroutine snoem(doy, kp, f107, z, mlat, nozm)

    real(sp) :: cosd, sind,thet
    sind(thet) = sin(thet/180.0*pi)
    cosd(thet) = cos(thet*pi/180.0)
    
    integer,intent(in)   :: doy
    real(sp),intent(in)  :: kp, f107
    real(sp),intent(out) :: z(16), mlat(33), nozm(33,16)

    real(sp) :: no_mean(33,16), eofs(33,16,3)
    real(sp) theta0         ! day number in degrees
    real(sp) dec            ! solar declination angle
    real(sp) m1, m2, m3     ! coefficients for first 3 eofs

    integer j, k, n


    !... read eof file

    open(unit=1,file='snoem_eof.dat',status='old')
        read(1,*) (z(k),k=1,16)
        read(1,*) (mlat(j),j=1,33)
        read(1,*) ((no_mean(j,k),j=1,33),k=1,16)
        read(1,*) (((eofs(j,k,n),j=1,33),k=1,16),n=1,3)
    close(unit=1)

    !... calculate coefficients (m1 to m3) for eofs based on geophysical parameters
    !... eof1 - kp 

    m1 =  kp * 0.689254 - 1.53366

    !... eof2 - declination

    theta0 = 360. * float(doy - 1) / 365.

    dec = 0.006918 &
        - 0.399912 * cosd(theta0)   + 0.070257 * sind(theta0)   &
        - 0.006758 * cosd(2*theta0) + 0.000907 * sind(2*theta0) &
        - 0.002697 * cosd(3*theta0) + 0.001480 * sind(3*theta0) 

    dec = dec * 180./3.1415927

    m2 = -0.31978  &
       + dec    * 0.097309   &
       + dec**2 * 0.00048979 &
       - dec**3 * 0.00010360 

    !... eof3 - f107 

    m3 =  alog10(f107) * 6.35777 - 13.8163 

    !... zonal mean distrib. is sum of mean and eofs

      nozm = no_mean                 &
                  - m1 * eofs(:,:,1) &
                  + m2 * eofs(:,:,2) &
                  - m3 * eofs(:,:,3) 

  end subroutine snoem
end module snoemmod
