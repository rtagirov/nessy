      module MOD_MSUB
      contains
      SUBROUTINE MSUB (A,B,JMAX,NP)
C***  A := A - B
      implicit real*8(a-h,o-z)

      DIMENSION A(NP,NP),B(NP,NP)
      DO 1 I=1,JMAX
      DO 1 K=1,JMAX
    1 A(I,K)=A(I,K)-B(I,K)
      RETURN
      end subroutine
      end module