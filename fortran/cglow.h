! This replaces glow.h:
!
! JMAX (number of vertical levels) is set to 64 for tgcm runs
! If not a tgcm run, JMAX is set in the namelist input file
!
!  integer :: JMAX
      integer,parameter :: lmax=123
      integer,parameter :: NMAJ=3
      integer,parameter :: NEI=10
      integer,parameter :: NEX=20
      integer,parameter :: NW=20
      integer,parameter :: NC=10
      integer,parameter :: NST=6
      integer,parameter :: NF=4

      integer, parameter :: dp = kind(1.0D0)
!      integer, parameter :: sp = kind(1.0)
      real(kind=dp), parameter    :: pi = 4.*ATAN(1.)
!radius of earth in centimeters
      real(kind=dp), parameter    :: Re = 6.371e8  
!gravitational constant dynes
      real(kind=dp), parameter    :: G  = 978.1    

! Standard parameters for photoelectron or aurora runs (up to 50 keV):
      integer, PARAMETER :: JMAX=120
      integer, PARAMETER :: NBINS=190
!
! Parameters for high energy aurora (up to 100 MeV):
!      integer,PARAMETER :: JMAX=170
!      integer,PARAMETER :: NBINS=343

