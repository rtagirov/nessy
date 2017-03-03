      module math

      implicit none

      private

      public :: expint
      public :: dsq, isq
      public :: dcube, icube
      public :: extrap_to_boundary

      contains

      RECURSIVE FUNCTION EXPINT(n, x) RESULT (Enx)

!     Calculates the exponential integral of order n

!     Using formulas 5.1.14, 5.1.53 and 5.1.56 from Handbook of Mathematical Functions,
!     Applied mathematics series 55 (1964) (ed. Abramovitz and Stegun)

      INTEGER, INTENT(IN) :: n

      REAL*8, INTENT(IN) :: x

      REAL*8 :: Enx

      IF (n .EQ. 0) THEN

         Enx = DEXP(-x) / x

      ELSEIF (n .EQ. 1) THEN

         Enx = E1(x) ! 5.1.53 and 5.1.56

      ELSE

         Enx = (DEXP(-x) - x * EXPINT(n - 1, x)) / (n - 1) ! 5.1.14

      ENDIF

      RETURN

      END FUNCTION EXPINT


      REAL*8 FUNCTION E1(x)

!     Calculates the exponential integral of order 1

      REAL*8, INTENT(IN) :: x

      REAL*8, DIMENSION(6) :: a

      REAL*8, DIMENSION(4) :: b, c

      DATA a /-0.57721566D0, 0.99999193D0, -0.24991055D0, 0.05519968D0, -0.00976004D0, 0.00107857D0/

      DATA b /8.5733287401D0, 18.0590169730D0, 8.6347608925D0, 0.2677737343D0/

      DATA c /9.5733223454D0, 25.6329561486D0, 21.0996530827D0, 3.9584969228D0/

      IF (x .GE. 0.0D0 .AND. x .LE. 1.0D0) THEN ! 5.1.53

         E1 = -DLOG(x) + a(1) + a(2) * x + a(3) * x**2.0D0 + 
     $         a(4) * x**3.0D0 + a(5) * x**4.0D0 + a(6) * x**5.0D0

      ELSEIF (x .GE. 1.0D0) THEN ! 5.1.56

             E1 = x**4.0D0 + b(1) * x**3.0D0 + b(2) * x**2.0D0 + b(3) * x + b(4)

             E1 = E1 / (x**4.0D0 + c(1) * x**3.0D0 + c(2) * x**2.0D0 + c(3) * x + c(4))

             E1 = E1 / x / DEXP(x)

      ENDIF

      RETURN

      END FUNCTION E1


      REAL*8 FUNCTION DSQ(x)

      REAL*8, INTENT(IN) :: x

      DSQ = x * x

      RETURN

      END FUNCTION DSQ


      REAL*8 FUNCTION ISQ(m)

      INTEGER, INTENT(IN) :: m

      ISQ = DBLE(m * m)

      RETURN

      END FUNCTION ISQ


      REAL*8 FUNCTION DCUBE(x)

      REAL*8, INTENT(IN) :: x

      DCUBE = x * x * x

      RETURN

      END FUNCTION DCUBE


      REAL*8 FUNCTION ICUBE(m)

      INTEGER, INTENT(IN) :: m

      ICUBE = DBLE(m * m * m)

      RETURN

      END FUNCTION ICUBE

      function extrap_to_boundary(n, x, y, bm2, bm1, b) result(yb)

!     linearly extrapolates the boundary value (yb) of the array y which is a 
!     function of x using the two preceding values (y(bm1), y(bm2))

      integer, intent(in) ::              n

      real*8, dimension(n), intent(in) :: x, y

      integer, intent(in) ::              b, bm1, bm2

      real*8 ::                           c1, c2, yb

      c1 = (y(bm2) - y(bm1)) / (x(bm2) - x(bm1))

      c2 = y(bm1) - c1 * x(bm1)

      yb = c1 * x(b) + c2

      return

      end function

      END MODULE MATH
