      module MOD_VSUB
      contains
C**********  MODULNAME: VSUB      ******* 24/03/87  22.14.04.******     7 KARTEN
      SUBROUTINE VSUB(A,B,N)
      implicit real*8(a-h,o-z)

      DIMENSION A(N),B(N)
      DO 1 I=1,N
    1 A(I)=A(I)-B(I)
      RETURN
      END subroutine
      end module