      MODULE MOD_GEOMESH

      CONTAINS

      SUBROUTINE GEOMESH(RADIUS, ENTOT, T, P, Z, RSTAR, ND, NP)

      USE MOD_PGRID
      USE MOD_RGRIDM

!     THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN RADIUS, P AND Z
!     P and Z mesh is needed for the ray-by-ray solution of the radiative transfer equation in spherical symmetry

      INTEGER, INTENT(OUT) ::  ND, NP
      REAL*8,  INTENT(IN)  ::  RSTAR

      real*8, allocatable, dimension(:), intent(out) :: Z, P
      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius

      CALL RGRIDM(RADIUS, ENTOT, T, RSTAR, ND)

      CALL PGRID(NP, ND, RADIUS, P)

      if (allocated(z)) deallocate(z)

      allocate(z(nd * np))

      z(1 : nd * np) = 0.0d0

      do L = 1, ND

         JMAX = NP + 1 - L

         do j = 1, JMAX

            i = (j - 1) * ND + L

            Z(i) = SQRT(radius(l)**2.0d0 - P(j)**2.0d0)

            write(*,'(A,2x,3(i4,2x),3(e15.7,2x))') 'geomesh:', l, j, i, radius(l)**2.0d0, P(j)**2.0d0, z(i)

        enddo

      enddo

      return

      end subroutine

      end module
