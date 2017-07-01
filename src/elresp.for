      module elresp ! ELectron RESPonse

      implicit none

      contains

      function ecd(L, ne, T) result(deriv)

      use vardatom_lte

      integer, intent(in) :: L ! depth index

      real*8,  intent(in) :: ne, T

      real*8              :: deriv, m1, m2, m3, m4

      integer             :: k, i, f

      deriv = 0.0d0

      if (natom_lte == 0) return

      do k = 1, natom_lte

         f  = 0
         m1 = 0.0
         m2 = 0.0
         m3 = 0.0
         m4 = 0.0

         do i = 1, lis_num

            if (lis_anum(i) /= eleatnum_lte(k)) cycle

            m1 = m1 + f * sef(i, ne, T)

            m2 = m2 + lis_cnum(i) * sef(i, ne, T)

            m3 = m3 + sef(i, ne, T)

            m4 = m4 + (1 + f) * lis_cnum(i) * sef(i, ne, T)

            f = f + 1

         enddo

         deriv = deriv + abxyzn_lte(k, L) * m4 * (m1 * m2 / m3 / m4 - 1.0d0) / m3

      enddo

      deriv = deriv / ne**2.0

      end function

      recursive function sef(i, ne, T) result(G)

      use vardatom_lte

      integer, intent(in) :: i     ! ionization stage number

      real*8, intent(in)  :: ne, T ! electron concentration and temperature

      real*8              :: G

      if (lis_cnum(i) == 0) then

         G = 1

         return

      else

         G = sef(i - 1, ne, T) * pfm(i, T) / pf(i - 1, T) / ne

      endif

      end function

      function pf(i, T) result(U)

      use phys
      use vardatom_lte

      integer, intent(in) :: i ! ionization stage number

      real*8,  intent(in) :: T ! temperature

      real*8              :: U

      integer             :: l

      U = 0

      do l = 1, lis_lnum(i)

         U = U + lis_weight(i, l) * dexp(-planck * light_speed * lis_levien(i, l) / boltz / T)

      enddo

      end function

      function pfm(i, T) result(V)

      use phys
      use vardatom_lte

      integer, intent(in) :: i ! ionization stage number

      real*8,  intent(in) :: T ! temperature

      real*8              :: V

      integer             :: l

      real*8              :: le

      le = planck / sqrt(2.0 * pai * elec_mass * boltz * T) ! de Broglie wavelength of electron

      V = 2.0 * pf(i, T) * dexp(-planck * light_speed * lis_levien(i - 1, 1) / boltz / T)

      V = V / le**3.0

      end function

      end module
