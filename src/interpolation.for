      module interpolation
      contains
      ! From Numerical Recipes Chapter 3 - Interpolation and Extrapolation - polint
      ! Polynomial Interpolation using Nevilles algorithm.
      !ms: f90tified it, +implicit none, +intent(*)
      !   Given arrays xa and ya, each of length n, and
      !   given a value x, this routine returns a value y, and
      !   an error estimate dy. If P (x) is the polynomial of
      !   degree N − 1 such that P (xai ) = yai ,
      !   i = 1, . . . , n, then the returned value y = P (x).
      SUBROUTINE polint(xa,ya,x,y,dy)
      implicit none
      !integer,intent(in) :: n
      real*8,intent(in)  :: x,xa(:),ya(:)
      real*8,intent(out) :: dy,y
      !integer,parameter  :: NMAX=10 ! Largest anticipated value of n.
      integer :: i,m,ns,n
      real*8 :: den,dif,dift,ho,hp,w,c(size(xa)),d(size(xa))
      ns=1
      n=size(xa)
      dif=abs(x-xa(1))
      do i=1,n                !  Here we ﬁnd the index ns of the closest table entry,
          dift=abs(x-xa(i))
          if (dift.lt.dif) then
                ns=i
                dif=dift
          endif
          c(i)=ya(i)             !  and initialize the tableau of c’s and d’s.
          d(i)=ya(i)
      enddo
      y=ya(ns)                    !  This is the initial approximation to y.
      ns=ns-1
      L13:do m=1,n-1              !  For each column of the tableau,
          do i=1,n-m         !  we loop over the current c’s and d’s and update them.
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                !  This error can occur only if two input xa’s are (to within roundoff) identical.
                if(den.eq.0.) stop 'polint: failure: xa not unique'
                den=w/den
                d(i)=hp*den           !  Here the c’s and d’s are updated.
                c(i)=ho*den
          enddo
          if (2*ns.lt.n-m)then   ! After each column in the tableau is completed, we decide
                dy=c(ns+1)        !  which correction, c or d, we want to add to our accu-
          else                   !  mulating value of y, i.e., which path to take through
                dy=d(ns)          !  the tableau—forking up or down. We do this in such a
                ns=ns-1           !  way as to take the most “straight line” route through the
          endif                  !  tableau to its apex, updating ns accordingly to keep track
          y=y+dy                 !  of where we are. This route keeps the partial approxima-
      enddo L13                   !  tions centered (insofar as possible) on the target x. The
      return                      !  last dy added is thus the error indication.
      END SUBROUTINE POLINT

      FUNCTION POLINT_F(xa,ya,x) result(y)
      real*8,intent(in) :: xa(:),ya(:),x
      real*8 :: y,dy
      call polint(xa,ya,x,y,dy)

      END FUNCTION
      end module

      !PROGRAM MAIN
      !   use interpolate
      !   integer,parameter :: n=3
      !   real*8 :: xa(n),ya(n),x,y,dy
      !   do i=1,n
      !     xa(i)=i;
      !   enddo
      !   ya(1)=1; ya(2)=4; ya(3)=9; ya(4)=16;
      !
      !   x=2.0; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=3.0; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=2.1; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=2.2; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=2.3; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=2.4; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !   x=2.5; call polint(xa,ya,n,x,y,dy);  print *,x,y,dy
      !end PROGRAM


