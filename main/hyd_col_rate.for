      MODULE HYD_COL_RATE

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: HYDCOLRATE

      CONTAINS

      FUNCTION HYDCOLRATE(nl, nu, T) RESULT(ColRate)

      USE MATH
      USE PHYS

!     e - H collisional rates according to L. C. Johnson, "Approximations For Collisional
!     And Radiative Transition Rates In Atomic Hydrogen", 1972, ApJ, 174, 227

!     The procedure works for ionization rates as well, in this case nu is equal to 11 and is not the principal quantum number

      INTEGER, INTENT(IN) :: nl, nu ! Principal numbers of upper and lower levels

      REAL*8, INTENT(IN) :: T ! Temperature

      REAL*8 :: ColRate ! Resulting collision rate (in cgs units)

      REAL*8 :: rnlnu

      REAL*8, DIMENSION(3) :: g

      REAL*8 :: gf

      REAL*8 :: Anlnu, Bnlnu

      REAL*8 :: Anl, Bnl

      REAL*8 :: C1, C2, C3, C4 ! Constant factors in the corresponding formulas

      REAL*8 :: B1, B2, B3 ! Bracket factors in the corresponding formulas

      REAL*8 :: x, y, z

      REAL*8 :: ynl, znl

      C1 = DSQRT(8.0D0 * boltz * T / pai / elec_mass)

      C4 = 32.0D0 / 3.0D0 / DSQRT(3.0D0) / pai

      g(1 : 3) = gfc(nl)

      IF (nu .NE. 11) THEN ! nl -> nu excitation collisional coefficient

         x = 1.0D0 - (ISQ(nl) / ISQ(nu))

         y = x * DABS(HYD_LEV_ENERGY(nl)) / boltz / T

         rnlnu = r(nl) * x

         z = rnlnu + y

         C2 = 2.0D0 * ISQ(nl) / x

         C3 = pai * DSQ(bohr_rad * y)

         B1 = (1.0D0 / y + 0.5D0) * EXPINT(1, y) - (1.0D0 / z + 0.5D0) * EXPINT(1, z)

         gf = g(1) + g(2) / x + g(3) / DSQ(x)

         Anlnu = C2 * C4 * nl * gf / DCUBE(nu * x)

         Bnlnu = (DSQ(C2) / ICUBE(nu)) * (1.0D0 + 4.0D0 / 3.0D0 / x + b(nl) / DSQ(x))

         B2 = Bnlnu - Anlnu * DLOG(C2)

         B3 = EXPINT(2, y) / y - EXPINT(2, z) / z

         ColRate = C1 * C2 * C3 * (Anlnu * B1 + B2 * B3)! * DEXP(y) / DSQRT(T)

      ELSE ! Ionization collisional coefficient

         ynl = DABS(HYD_LEV_ENERGY(nl)) / boltz / T

         znl = r(nl) + ynl

         C2 = 2.0D0 * pai * DSQ(nl * bohr_rad * ynl)

         B1 = EXPINT(1, ynl) / ynl - EXPINT(1, znl) / znl

         Bnl = 2.0D0 * ISQ(nl) * (5.0D0 + b(nl)) / 3.0D0

         Anl = C4 * nl * (g(1) / 3.0D0 + g(2) / 4.0D0 + g(3) / 5.0D0)

         B2 = Bnl - Anl * DLOG(2.0D0 * ISQ(nl))

         B3 = ksi(ynl) - ksi(znl)

         ColRate = C1 * C2 * (Anl * B1 + B2 * B3)! * DEXP(ynl) / DSQRT(T)

      ENDIF

      RETURN

      END FUNCTION HYDCOLRATE

 
      REAL*8 FUNCTION r(n)

      INTEGER, INTENT(IN) :: n

      IF (n .EQ. 1) THEN; r = 0.45D0; ELSE; r = 1.94D0 * DBLE(n)**(-1.57D0); ENDIF

      RETURN

      END FUNCTION r


      FUNCTION gfc(n) RESULT(g)

      INTEGER, INTENT(IN) :: n

      REAL*8, DIMENSION(3) :: g

      IF (n .EQ. 1) THEN

         g = (/1.1330D0, -0.4059D0, 0.07014D0/)

      ELSEIF (n .EQ. 2) THEN

         g = (/1.0785D0, -0.2319D0, 0.02947D0/)

      ELSE

         g(1) = 0.9935D0 + 0.2328D0 / n - 0.1296D0 / n / n

         g(2) = (-0.6282D0 + 0.5598D0 / n - 0.5299D0 / n / n) / n

         g(3) = (0.3887D0 - 1.181D0 / n + 1.470D0 / n / n) / n / n

      ENDIF

      RETURN

      END FUNCTION gfc


      REAL*8 FUNCTION b(n)

      INTEGER, INTENT(IN) :: n

      IF (n .EQ. 1) THEN

          b = -0.603D0

      ELSE

          b = (4.0D0 - 18.63D0 / n + 36.24D0 / n / n - 28.09D0 / n / n / n) / n

      ENDIF

      RETURN

      END FUNCTION b


      REAL*8 FUNCTION ksi(t)

      USE MATH

      REAL*8, INTENT(IN) :: t

      ksi = EXPINT(0, t) - 2.0D0 * EXPINT(1, t) + EXPINT(2, t)

      RETURN

      END FUNCTION ksi

      END MODULE HYD_COL_RATE
