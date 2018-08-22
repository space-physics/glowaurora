      subroutine vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1,lat)
      implicit none
!
! Determines five roots of the equation:
!   a4*x**4 + 4.*a3*x**3 + 6.*a2*x**2 + 4.*a1*x + a0 = 0.
!
! Procedure is specificlly designed for real quartics with real roots
! only one of which is positive.
!
! This is called by elden for electron density.
! 
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in)  :: a0,a1,a2,a3,a4
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: root
!
! Local:
      integer :: k,i,nlevs,i0,i1
      real,dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3 ! work arrays
      real,parameter :: e=1.e-300 ! largest exponent on ieee is about 307
!
      nlevs = lev1-lev0+1
      i0 = lon0 ; i1 = lon1
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! w1 = ch
          w1(k,i) = -(a4(k,i)*a0(k,i)-4.*a3(k,i)*a1(k,i)+3.*a2(k,i)**2)/
     |                12.
!
! w2 = cg
          w2(k,i) = (a4(k,i)*(a2(k,i)*a0(k,i)-a1(k,i)**2)-a3(k,i)*
     |      (a3(k,i)*a0(k,i)-a1(k,i)*a2(k,i))+a2(k,i)*(a3(k,i)*a1(k,i)-
     |      a2(k,i)**2))/4.
!
! root=rlam=-2.*real((.5*(cmplx(cg,0.)+csqrt(cmplx(cg**2+4.
!      *ch**3+e,0.)))+cmplx(e,0.))**(1./3.))
! 
          root(k,i) = -2.*real((.5*(cmplx(w2(k,i),0.)+
     |      csqrt(cmplx(w2(k,i)**2+4.*w1(k,i)**3+e,0.)))+
     |      cmplx(e,0.))**(1./3.))
!
! W1=P=SQRT(A(5)*RLAM+A(4)**2-A(5)*A(3)+E)
!
          w1(k,i) = a4(k,i)*root(k,i)+a3(k,i)**2-a4(k,i)*a2(k,i)+e
          if (w1(k,i) < 0.) w1(k,i) = 0.
          w1(k,i) = sqrt(w1(k,i))
!
! W2=Q=SQRT((2.*RLAM+A(3))**2-A(5)*A(1)+E)
!
          w2(k,i) = sqrt((2.*root(k,i)+a2(k,i))**2-a4(k,i)*a0(k,i)+e)
!
! W3=PQ=2.*A(4)*RLAM+A(4)*A(3)-A(5)*A(2)+E
!
          w3(k,i) = 2.*a3(k,i)*root(k,i)+a3(k,i)*a2(k,i)-a4(k,i)*a1(k,i)
     |      +e
!
!  W1=P=SIGN(P,Q*PQ)
!
          w1(k,i) = sign(w1(k,i),w2(k,i)*w3(k,i))
!
! W3=P-A4
!
          w3(k,i) = w1(k,i)-a3(k,i)
!
! Final evaluation of root:
!
          root(k,i) = (w3(k,i)+sqrt(w3(k,i)**2-a4(k,i)*(a2(k,i)+2.*
     |      root(k,i)-w2(k,i))))/a4(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1


      end subroutine vquart
