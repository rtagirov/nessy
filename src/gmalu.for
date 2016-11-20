      module MOD_GMALU
      contains
C**********  MODULNAME: GMALU     ******* 24/03/87  19.29.37.******     9 KARTEN
      SUBROUTINE GMALU (GA,U,V,LMAX)
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION GA (LMAX),U(LMAX),V(LMAX)
      LZ=LMAX-1
      DO 1 L=1,LZ
    1 V(L)=GA(L)*(U(L)-U(L+1))
      RETURN
      END subroutine
      end module
