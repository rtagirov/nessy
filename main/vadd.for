      module MOD_VADD
      contains
      SUBROUTINE VADD (A,B,N)
C***  VECTOR ADDITION  A = A + B
      implicit real*8(a-h,o-z)

      DIMENSION A(N),B(N)
      DO 1 I=1,N
    1 A(I)=A(I)+B(I)
      RETURN
      end subroutine
      end module