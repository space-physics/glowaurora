! GLOW header file sets altitude, electron energy, and solar grid sizes.
!
! Comment out all but one of the following before compiling:
!
! Standard parameters for photoelectron or aurora runs (up to 50 keV):
      integer, PARAMETER :: JMAX=120
      integer, PARAMETER :: NBINS=190
      integer, PARAMETER :: LMAX=123
!
! Parameters for high energy aurora (up to 100 MeV):
!     PARAMETER (JMAX=170)
!     PARAMETER (NBINS=343)
!     PARAMETER (LMAX=123)
