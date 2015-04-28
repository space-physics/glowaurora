MODULE machprec
    public
      integer, parameter :: dp = kind(1.0D0)
      integer, parameter :: sp = kind(1.0)
      real(sp), parameter:: pi = 4*ATAN(1.)
END MODULE machprec
