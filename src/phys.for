      MODULE PHYS

      IMPLICIT NONE

      REAL*8, PARAMETER :: pai = 3.1415926

      REAL*8, PARAMETER :: boltz = 1.3806488D-16       ! Boltzmann constant, erg / K

      REAL*8, PARAMETER :: elec_mass = 9.1093829D-28   ! Electron mass, g

      REAL*8, PARAMETER :: prot_mass = 1.6726217D-24   ! Proton mass, g

      REAL*8, PARAMETER :: bohr_rad = 5.29D-9          ! Bohr radius, cm

      REAL*8, PARAMETER :: fine_str = 7.2973525D-3     ! Fine structure constant

      REAL*8, PARAMETER :: light_speed = 2.9979245D+10 ! Speed of light, cm / s

      REAL*8, PARAMETER :: planck = 6.6260688D-27      ! Planck constant, erg * s

      real*8, parameter :: elec_charg = 4.80320425d-10 ! Electron charge in statcoulombs

      real*8, parameter :: sigma_thomson = 6.65d-25    ! Thomson cross-section in cm^{-2}

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

      pure function BNUE(XLAMBDA, T)

!     PLANCK FUNCTION, LAMBDA IN ANGSTROEM, T IN KELVIN
!     BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ
!     CONSTANTS : C1 = H * C / K (DIMENSION ANGSTROEM * KELVIN )
!                 C2 = 2 * H * C (DIMENSION ANGSTROEM**3 * ERG/SEC/HZ/CM**2)

      implicit none

      real*8 :: BNUE, cxt
      real*8, intent(in) :: XLAMBDA, T
      real*8, parameter  :: C1 = 1.438831D8
      real*8, parameter  :: C2 = 3.972967D+8
  
      if (T .le. 0.0d0) then

         BNUE = 0.0d0

      else

         CXT = C1 / XLAMBDA / T

        if (CXT .gt. 500.0d0) then

           BNUE = 0.0

        else

           BNUE = C2 / (EXP(CXT) - 1d0) / XLAMBDA**3

        endif

      endif

      end function

      function opt_dep(opac, h, ND) result(tau)

!     Calculate optical depth (tau) from a given 
!     opacity (opac) on a given height grid (h)

      integer, intent(in) :: ND

      real*8, intent(in), dimension(ND) :: opac, h

      real*8, dimension(ND) :: tau

      integer :: l

      tau(1) = 0.0D0

      do l = 2, ND

         tau(l) = (h(l - 1) - h(l)) * (opac(l - 1) + opac(l)) / 2.0D0

      enddo

      return

      end function opt_dep

      function sigma_rayleigh(wvl) result(sigma)

      use vardatom

      real*8, intent(in) :: wvl ! wavelength in angstroems

      real*8 :: sigma

      integer :: u, l

      real*8 :: f, d, r, w

      real*8 :: wcm, wA

      integer :: j, n_hyd_lev, n_hyd_lin

      n_hyd_lev = nlast(1) - 2

      n_hyd_lin = 0

      do j = 1, n_hyd_lev - 1; n_hyd_lin = n_hyd_lin + j; enddo

!      print*, 'n_hyd_lin = ', n_hyd_lin

      sigma = 0.0d0

      do j = 1, n_hyd_lin

         u = indnup(j); l = indlow(j)

         if (l .eq. 2) then

             wcm = 1.0d0 / (elevel(u) - elevel(l))

             wA = 1.0d8 * wcm

             if (abs(wvl - wA) .lt. 1.0d0) cycle

             f = (elec_mass * light_speed * wcm**2.0d0 * weight(l) * einst(u, l)) /
     $           (8.0d0 * pai**2.0d0 * elec_charg**2.0d0 * weight(u))

             r = wvl / wA

             d = (r - 1.0d0) * (r + 1.0d0)

             sigma = sigma + f / d

         endif

      enddo

      sigma = sigma_thomson * sigma**2.0d0

      return

      end function

      END MODULE PHYS
