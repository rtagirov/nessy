      MODULE MOD_STORXJC

      CONTAINS

      SUBROUTINE STORXJC(XJCREA, XJC, EDDREA, EDDI, nd, nf, K)
!     store XJC and EDDI for the frequency K

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION xjcrea(nd, nf), xjc(nd), eddrea(3, nd, nf), eddi(3, nd)

      DO L = 1, nd

         xjcrea(L, K) =    xjc(L)
         eddrea(1, L, K) = eddi(1, L)
         eddrea(2, L, K) = eddi(2, L)
         eddrea(3, L, K) = eddi(3, L)

      ENDDO

      RETURN

      END SUBROUTINE

      END MODULE
