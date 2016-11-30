      MODULE MOD_RGRIDM

      CONTAINS

      SUBROUTINE RGRIDM(radius, entot, T, rstar, ND)

      USE MOD_ERROR
      USE COMMON_BLOCK

      use file_operations

      real*8,  intent(in)  :: rstar

      integer, intent(out) :: ND

      real*8, allocatable, dimension(:), intent(out) :: T, entot, radius

      real*8 :: vt, elect, v1, v2, v4

!     T:      Temperature in K
!     elect:  ELECTRON DENSITY IN ELECTRON / cm^3
!     entot:  HEAVY PARTICLE DENSITY
!     vt:     Turbulent velocity in km / s
!     radius: HEIGHT IN UNITS OF SOLAR RADII

      ND = num_of_lines('FAL_VD'); DPN = ND

      if (allocated(T))      deallocate(T)
      if (allocated(entot))  deallocate(entot)
      if (allocated(radius)) deallocate(radius)
      if (allocated(height)) deallocate(height)

      allocate(T(ND))
      allocate(entot(ND))
      allocate(radius(ND))
      allocate(height(ND))

      open(9, file = 'FAL_VD')

      do i = 1, ND

         read(9, *) v1, v2, elect, v4, vt

         height(i) = v1
         T(i) =      v2
         entot(i) =  v4

         radius(i) = 1.0d0 + height(i) * 1.0d5 / rstar ! height in km, rstar and radius in cm

      enddo

      close(9)

      end subroutine

      end module
