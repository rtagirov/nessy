      MODULE MOD_RGRIDM

      CONTAINS

      SUBROUTINE RGRIDM(RADIUS, ENTOT, T, RMAX, RSTAR, AMU, ATMEAN, ND)

      USE MOD_ERROR
      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      integer,     intent(out) :: ND

      real*8,      intent(in)  :: RMAX, RSTAR, AMU, ATMEAN

      real*8, allocatable, dimension(:), intent(out) :: T, ENTOT, RADIUS

      DIMENSION T(NDDIM),ENTOT(NDDIM),RADIUS(NDDIM),height(nddim),VDTAB(NDDIM)

!     ENTOT: HEAVY PARTICLE DENSITY
!     XNETAB(L): ELECTRON DENSITY IN ELECTRON/CM^3
!     RSUN = 6.960E^10 CM
!     RADIUS = HEIGHT IN UNITS OF SOLAR RADII

      ND = 1

      VDTAB(1:NDDIM) = 0

      OPEN(9, FILE = 'FAL_VD', STATUS = 'OLD')

  400 READ(9, *, end = 55) height(ND), T(ND), XNETAB(ND), ENTOT(ND), VDTAB(ND)

              RADIUS(ND) = 1.D0 + height(ND) * 1.D5 / RSTAR  ! HEIGTH in km, RSTAR in cm

          ND = ND+1
          GOTO 400

   55     CLOSE (9)
          ND=ND-1
          CONTINUE

      DPN = ND

      END subroutine

      end module
