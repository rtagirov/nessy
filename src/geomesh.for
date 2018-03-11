      module geo_mesh

      implicit none

      contains

      subroutine geomesh(radius, entot, T, p, z, rstar, amu, atmean, ND, NP)

      use common_block
      use file_operations

!     THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN RADIUS, P AND Z
!     P and Z mesh is needed for the ray-by-ray solution of the radiative transfer equation in spherical symmetry

      integer, intent(out) ::  ND, NP
      real*8,  intent(in)  ::  rstar

      real*8, intent(in) :: amu, atmean

      real*8, allocatable, dimension(:), intent(out) :: Z, P
      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius

      logical :: atm_mod_file_exists

      inquire(file = atm_mod_file, exist = atm_mod_file_exists)

      if (.not. atm_mod_file_exists) stop 'Atmosphere model file has not been found. Abort.'

!     height: height in km
!     T:      Temperature in K
!     entot:  heavy particle density
!     radius: height in units of solar radii

      selectcase(num_of_columns(atm_mod_file))

          case(4);  call read_mur_mod(rstar, height, radius, T, entot, ND)              ! Atmosphere model in MURAM format

          case(5);  call read_fal_mod(rstar, height, radius, T, entot, ND)              ! Atmosphere model in FAL format

          case(7);  call read_kur_mod(amu, atmean, rstar, height, radius, T, entot, ND) ! Atmosphere model in Kurucz 7  columns format

          case(10); call read_kur_mod(amu, atmean, rstar, height, radius, T, entot, ND) ! Atmosphere model in Kurucz 10 columns format

          case default; stop 'Function read_atm_file_col: col argument is not recognized. Abort.'

      endselect

      NDPMIN = tempmin(T, ND)

      call pgrid(NP, ND, radius, p)

      call zgrid(radius, p, z, ND, NP)

      call print_strat(ND, height, T, entot, radius)

      dpn = ND

      return

      end subroutine

      subroutine print_strat(ND, h, T, n, r)

      use math

      integer, intent(in) :: ND

      real*8, dimension(ND), intent(in) :: h, T, n, r

      real*8, dimension(ND) :: dT, dh, dn, dr, gradT, gradn

      integer :: fu, i

      character(len = 1000) :: fmt_head, fmt_body

      fmt_head = '(A,8x,A,13x,A,11x,A,14x,A,13x,A,14x,A,8x,A,/)'

      fmt_body = '(i3,4x,2(F9.2,4x),es9.3,4x,F9.2,4x,es15.7,4x,F9.2,4x,es10.3))'

      do i = 1, ND - 1

         dh(i) = h(i + 1) - h(i)

         dr(i) = r(i + 1) - r(i)

         dT(i) = T(i + 1) - T(i)

         dn(i) = n(i + 1) - n(i)

      enddo

      dh(ND) = extrap_to_boundary(ND, h, dh, 1)
      dr(ND) = extrap_to_boundary(ND, h, dr, 1)
      dT(ND) = extrap_to_boundary(ND, h, dT, 1)
      dn(ND) = extrap_to_boundary(ND, h, dn, 1)

      gradT = dT / dh

      gradn = dn * maxval(abs(dr)) / dr / maxval(abs(dn))

      fu = 1242; open(unit = fu, file = 'strat.out', action = 'write')

      write(fu, fmt_head) 'idx', 'h', 'T', 'n', 'dh', 'r', 'dT/dh', 'dn/dr'

      do i = 1, ND; write(fu, fmt_body) i, h(i), T(i), n(i), abs(dh(i)), r(i), gradT(i), gradn(i); enddo

      close(fu)

      end subroutine

      SUBROUTINE PGRID(NP, ND, R, P)
!     GRID OF IMPACT-PARAMETER POINTS

      integer,                intent(in) ::  ND
      real*8,  dimension(ND), intent(in) ::  R

      integer, intent(out) :: NP

      real*8, allocatable, dimension(:), intent(out) :: P

      real*8, dimension(13) :: cp

      integer :: NC, L, J

      data cp /0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 
     $         0.91652d0, 0.97980d0, 0.99499d0, 0.99875d0/

!     NC = NUMBER OF CORE-INTERSECTING RAYS

      NC = 13; NP = ND + NC

      allocate(P(NP))

      P(1 : NC) = cp(1 : NC)

!     CORE RAYS EQUELLY SPACED

      do L = 1, ND; J = NP + 1 - L; P(J) = R(L); enddo

      return

      end subroutine

      SUBROUTINE ZGRID(RADIUS, P, Z, ND, NP)

!     THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN Z

      integer, intent(in)  :: ND, NP
      real*8,  intent(in)  :: RADIUS(ND), P(NP)

      real*8,  allocatable, dimension(:), intent(out) :: Z

      integer :: i, j, L, JMAX

      allocate(z(nd * np))

      z(1 : nd * np) = 0.0d0

      do L = 1, ND

         JMAX = NP + 1 - L

         do j = 1, JMAX

            i = (j - 1) * ND + L

            Z(i) = SQRT((radius(l) - P(j)) * (radius(l) + P(j)))

        enddo

      enddo

      return

      end subroutine

      subroutine read_mur_mod(rstar, h, r, T, n, ND)

      use phys
      use file_operations

      real*8,  intent(in)  :: rstar

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, n, r, h

      real*8, allocatable, dimension(:) :: p

      ND = num_of_lines(atm_mod_file)

      allocate(h(ND))
      allocate(T(ND))
      allocate(n(ND))
      allocate(r(ND))
      allocate(p(ND))

      h = read_atm_file_col(1) ! Height
      T = read_atm_file_col(2) ! Temperature
      p = read_atm_file_col(3) ! Pressure

      n = p / (boltz * T) ! Number density

      h = abs(h - maxval(h))

      r = 1.0d0 + h * 1.0d5 / rstar ! Calculation of radius in relative units: h (height) in km, rstar in cm

      deallocate(p)

      end subroutine

      subroutine read_fal_mod(rstar, h, r, T, n, ND)

      use file_operations

      real*8,  intent(in)  :: rstar

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, n, r, h

      ND = num_of_lines(atm_mod_file)

      allocate(h(ND))
      allocate(T(ND))
      allocate(n(ND))
      allocate(r(ND))

      h = read_atm_file_col(1) ! Height
      T = read_atm_file_col(2) ! Temperature
      n = read_atm_file_col(4) ! Number density

      r = 1.0d0 + h * 1.0d5 / rstar ! Calculation of radius in relative units: h (height) in km, rstar in cm

      end subroutine

      subroutine read_kur_mod(amu, atmean, rstar, height, radius, T, entot, ND)

      use file_operations

      real*8,  intent(in)  :: rstar

      real*8,  intent(in)  :: amu, atmean

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius, height

      real*8, allocatable, dimension(:) :: entotn, delr

      real*8, allocatable, dimension(:) :: elec_conc, rho, vturb, pressure

      integer :: k, l

!     ak - Boltzmann constant
      real*8, parameter :: ak =  1.38062259d-16
      real*8, parameter :: mun = 1.66054d-24

      ND = num_of_lines(atm_mod_file)

      allocate(T(ND))
      allocate(entot(ND))
      allocate(radius(ND))
      allocate(height(ND))

      allocate(rho(ND))
      allocate(pressure(ND))
      allocate(elec_conc(ND))
      allocate(vturb(ND))

      allocate(entotn(ND))
      allocate(delr(ND - 1))

      rho =       read_atm_file_col(1)
      T =         read_atm_file_col(2)
      pressure =  read_atm_file_col(3)
      elec_conc = read_atm_file_col(4)
      vturb =     read_atm_file_col(7)

!     taking into account turbulent pressure
!      entotn = pressure / (ak * T + 0.5d0 * atmean * mun * vturb**2.0d0)
      entotn = pressure / ak / T

      entot  = entotn - elec_conc

      do l = 1, ND - 1

         delr(l) = (2 / (amu * atmean * rstar)) * (rho(l + 1) - rho(l)) / (entot(l + 1) + entot(l))

      enddo

      radius(ND) = 1.0

      do k = 1, ND - 1; radius(ND - k) = radius(ND - k + 1) + delr(ND - k); enddo

      height = (radius - radius(ND)) * rstar / 1D+5 ! height in km, rstar in cm

      deallocate(rho)
      deallocate(vturb)
      deallocate(elec_conc)
      deallocate(pressure)
      deallocate(delr)
      deallocate(entotn)

      end subroutine

      function tempmin(T, ND) result(NDPMIN)

      integer :: NDPMIN

      integer, intent(in) :: ND

      real*8, dimension(ND), intent(in) :: T(ND)

      integer :: L

      NDPMIN = 1

      do L = 2, ND - 1

         if ((T(L) .lt. T(L - 1)) .AND. (T(L) .lt. T(L + 1)) .and. T(L) .lt. T(NDPMIN)) NDPMIN = L

      enddo

      return

      end function

      end module
