      module MOD_MVMD
      contains
      SUBROUTINE MVMD (BX,B,C,JMAX,JMM,NP)
C***  MATRIX (VOLL)  B  *  MATRIX (DIAGONAL)  C
C***  ERGEBNIS-MATRIX  BX(VOLL)
C*** AKTUELLES FORMAT BX(JMAX,JMM)=B(JMAX,JMAX)*C(JMAX,JMM)
C***  WOBEI DIE UEBERZAEHLIGEN ZEILEN DER DIAGONALMATRIX C VERSCHWINDEN
      implicit real*8(a-h,o-z)

      DIMENSION BX(NP,NP),B(NP,NP),C(NP)
      DO 1 K=1,JMM
      CK=C(K)
      DO 1 I=1,JMAX
    1 BX(I,K)=B(I,K)*CK
      RETURN
      end subroutine
      end module