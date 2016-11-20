      module MOD_MDV
      contains
      SUBROUTINE MDV (A,W,N)
C*** MATRIX A (DIAGONAL)  *  VEKTOR W
C***  ERGEBNIS-VEKTOR UEBERSCHREIBT  W
      implicit real*8(a-h,o-z)

      DIMENSION  A(N),W(N)
      DO 1 I=1,N
    1 W(I)=A(I)*W(I)
      RETURN
      end subroutine
      end module