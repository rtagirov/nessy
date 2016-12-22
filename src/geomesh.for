      MODULE MOD_GEOMESH

      CONTAINS

      SUBROUTINE GEOMESH(RADIUS, ENTOT, T, P, Z, RSTAR, ND, NP)

!     THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN RADIUS, P AND Z
!     P and Z mesh is needed for the ray-by-ray solution of the radiative transfer equation in spherical symmetry

      INTEGER, INTENT(OUT) ::  ND, NP
      REAL*8,  INTENT(IN)  ::  RSTAR

      real*8, allocatable, dimension(:), intent(out) :: Z, P
      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius

      CALL RGRIDM(RADIUS, ENTOT, T, RSTAR, ND)

      CALL PGRID(NP, ND, RADIUS, P)

      CALL ZGRID(RADIUS, P, Z, ND, NP)

      return

      end subroutine

      SUBROUTINE RGRIDM(radius, entot, T, rstar, ND)

      USE MOD_ERROR
      USE COMMON_BLOCK
      use FILE_OPERATIONS

      real*8,  intent(in)  :: rstar

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius

      real*8 :: vt, elect, v1, v2, v4

!     height: height in km
!     T:      Temperature in K
!     entot:  HEAVY PARTICLE DENSITY
!     radius: HEIGHT IN UNITS OF SOLAR RADII

      ND = num_of_lines(atm_mod_file); DPN = ND

      if (allocated(T))      deallocate(T)
      if (allocated(entot))  deallocate(entot)
      if (allocated(radius)) deallocate(radius)
      if (allocated(height)) deallocate(height)

      allocate(T(ND))
      allocate(entot(ND))
      allocate(radius(ND))
      allocate(height(ND))

      height = read_atm_mod(atm_mod_file, '1')
      T =      read_atm_mod(atm_mod_file, '2')
      entot =  read_atm_mod(atm_mod_file, '4')

      radius = 1.0d0 + height * 1.0d5 / rstar ! height in km, rstar and radius in cm

      end subroutine

      SUBROUTINE PGRID(NP, ND, R, P)
!     GRID OF IMPACT-PARAMETER POINTS

      integer,                intent(in) ::  ND
      real*8,  dimension(ND), intent(in) ::  R

      integer, intent(out) :: NP

      real*8, allocatable, dimension(:), intent(out) :: P

      real*8, dimension(13) :: cp

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

      RETURN

      END subroutine

      end module
