      MODULE PHYS

      IMPLICIT NONE

      REAL*8, PARAMETER :: pai = 3.1415926

      REAL*8, PARAMETER :: boltz = 1.3806488D-16 ! Boltzmann constant, erg / K

      REAL*8, PARAMETER :: elec_mass = 9.1093829D-28 ! Electron mass, g

      REAL*8, PARAMETER :: prot_mass = 1.6726217D-24 ! Proton mass, g

      REAL*8, PARAMETER :: bohr_rad = 5.29D-9 ! Bohr radius, cm

      REAL*8, PARAMETER :: fine_str = 7.2973525D-3 ! Fine structure constant

      REAL*8, PARAMETER :: light_speed = 2.9979245D+10 ! Speed of light, cm / s

      REAL*8, PARAMETER :: planck = 6.6260688D-27 ! Planck constant, erg * s

      PUBLIC

      CONTAINS

      REAL*8 FUNCTION HYD_LEV_ENERGY(n) ! Energy of a hydrogen level with principal quantum number n (in CGS units)

      USE MATH

      INTEGER, INTENT(IN) :: n

      HYD_LEV_ENERGY = -elec_mass * DSQ(light_speed * fine_str)  / ISQ(n) / 2.0D0

      RETURN

      END FUNCTION HYD_LEV_ENERGY

      REAL*8 FUNCTION PLANCK_FUNC(nu, T) ! Planck intensity per unit frequency (in CGS units)

      USE MATH

      REAL*8, INTENT(IN) :: nu, T

      REAL*8 :: PLANCK_EXP

      REAL*8 :: PLANCK_FREQ_FACT

      PLANCK_FREQ_FACT = 2.0D0 * planck * DCUBE(nu) / DSQ(light_speed)

      PLANCK_EXP = DEXP(planck * nu / boltz / T)

      PLANCK_FUNC =  PLANCK_FREQ_FACT / (PLANCK_EXP - 1.0D0)

      RETURN

      END FUNCTION PLANCK_FUNC

      function optical_depth(opacity, height, ND) result(tau)

      integer, intent(in) :: ND

      real*8, intent(in), dimension(ND) :: opacity, height

      real*8, dimension(ND) :: tau

      integer :: l

      tau(1) = 0.0D0

      do l = 2, ND

         tau(l) = (height(l - 1) - height(l)) * (opacity(l - 1) + opacity(l)) / 2.0D0

      enddo

      return

      end function optical_depth

      END MODULE PHYS
