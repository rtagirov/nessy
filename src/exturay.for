      module MOD_EXTURAY

      contains

      subroutine extUray (U, Uray, nd, np, iray)
!     now store Uray for the p-ray
      IMPLICIT REAL*8(A - H, O - Z)

      dimension U(nd, np), Uray(nd)

      do L = 1, nd; Uray(L) = U(L, iray); enddo

      return

      end subroutine

      end module
