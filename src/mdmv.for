      module MOD_MDMV
      contains
      SUBROUTINE MDMV (A,B,JMAX,NP)
C***  MATRIX (DIAGONAL)  A  *  MATRIX (VOLL)  B
C***  ERGEBNIS-MATRIX UEBERSCHREIBT  B
      implicit real*8(a-h,o-z)

      DIMENSION A(NP),B(NP,NP)
      DO 1 I=1,JMAX
      AI=A(I)
      DO 1 K=1,JMAX
    1 B(I,K)=B(I,K)*AI
      RETURN
      end subroutine
      end module