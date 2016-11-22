      module MOD_ACOPY

      contains

      SUBROUTINE ACOPY (A1, A2, N)

      implicit real*8(a - h, o - z)

      dimension A1(N, N), A2(N, N)

      do J = 1, N

         do I = 1, N

            A1(J, I) = A2(J, I)

         enddo

      enddo

      return

      end subroutine

      end module
