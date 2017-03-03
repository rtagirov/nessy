      module geo_mesh

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

      integer :: fu, i

      inquire(file = atm_mod_file, exist = atm_mod_file_exists)

      if (.not. atm_mod_file_exists) stop 'Atmosphere model file has not been found. Abort.'

!     height: height in km
!     T:      Temperature in K
!     entot:  HEAVY PARTICLE DENSITY
!     radius: HEIGHT IN UNITS OF SOLAR RADII

      selectcase(num_of_columns(atm_mod_file))

          case(5); call read_fal_fmt_mod(rstar, height, radius, T, entot, ND)

          case(7); call read_kur_fmt_mod(amu, atmean, rstar, height, radius, T, entot, ND)

          case default; stop 'Function read_atm_file_col: col argument is not recognized. Abort.'

      endselect

      call pgrid(NP, ND, radius, p)

      call zgrid(radius, p, z, ND, NP)

      fu = 1242; open(unit = fu, file = 'ATM_STR', action = 'write')

      do i = 1, ND; write(fu, '(i3,4x,2(F9.2,4x),es15.7)') i, height(i), T(i), entot(i); enddo

      close(fu)

      dpn = ND

      return

      end subroutine

      SUBROUTINE PGRID(NP, ND, R, P)
!     GRID OF IMPACT-PARAMETER POINTS

      implicit none

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

      implicit none
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

      subroutine read_fal_fmt_mod(rstar, height, radius, T, entot, ND)

      use file_operations

      implicit none

      real*8,  intent(in)  :: rstar

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius, height

      ND = num_of_lines(atm_mod_file)

      if (allocated(height)) deallocate(height)

      allocate(T(ND))
      allocate(entot(ND))
      allocate(radius(ND))
      allocate(height(ND))

      height = read_atm_file_col(1)
      T =      read_atm_file_col(2)
      entot =  read_atm_file_col(4)

      radius = 1.0d0 + height * 1.0d5 / rstar ! height in km, rstar in cm

      end subroutine

      subroutine read_kur_fmt_mod(amu, atmean, rstar, height, radius, T, entot, ND)

      use file_operations

      implicit none

      real*8,  intent(in)  :: rstar

      real*8, intent(in) :: amu, atmean

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius, height

      real*8, allocatable, dimension(:) :: entotn, delr

      real*8, allocatable, dimension(:) :: elec_conc, rho, vturb, pressure

!     ak - Boltzmann constant
      real*8, parameter :: ak =  1.38062259d-16
      real*8, parameter :: MUN = 1.66054d-24

      ND = num_of_lines(atm_mod_file)

      if(allocated(height)) deallocate(height)

      allocate(height(ND))
      allocate(radius(ND))
      allocate(T(ND))
      allocate(entot(ND))

      allocate(rho(ND))
      allocate(pressure(ND))
      allocate(elec_conc(ND))
      allocate(vturb(ND))
      allocate(entotn(ND))
      allocate(delr(ND))

      rho =       read_atm_file_col(1)
      T =         read_atm_file_col(2)
      pressure =  read_atm_file_col(3)
      elec_conc = read_atm_file_col(4)
      vturb =     read_atm_file_col(7)

      pressure = rho * 10.0**4.44 ! this is from the old version and I do not understand what's it doing there
!      pressure = rho * 10.0**4.5 ! this is from the old version and I do not understand what's it doing there

!     TAKING INTO ACCOUNT TURBULEN PRESSURE
              
      entotn = pressure / (AK * T + 0.5 * ATMEAN * MUN * vturb**2.)

      entot  = entotn - elec_conc

      do l = 1, ND - 1

         delr(l) = (2 / (amu * atmean * rstar)) * (rho(l + 1) - rho(l)) / (entot(l + 1) + entot(l))

      enddo

      radius(ND) = 1.0

      do k = 1, ND - 1; radius(ND - k) = radius(ND - k + 1) + delr(ND - k); enddo

      height = (radius - radius(ND)) * rstar / 1D+5 ! height in km, rstar in cm

      deallocate(rho)
      deallocate(pressure)
      deallocate(elec_conc)
      deallocate(vturb)
      deallocate(entotn)
      deallocate(delr)

      end subroutine

      end module
