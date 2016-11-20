      module MOD_MVV
      contains
      SUBROUTINE MVV (WX,B,W,JMAX,JMM,NP)
C***  MATRIX (VOLL)  B  *  VEKTOR W
C***  ERGEBNIS-VEKTOR  WX
C***  AKTUELLES FORMAT  WX(JMAX) = B(JMAX,JMM) * W(JMM)
      implicit real*8(a-h,o-z)

      DIMENSION WX(NP),B(NP,NP),W(NP)
      DO 1 I=1,JMAX
      WXI = .0
      DO 2 K=1,JMM
    2 WXI=WXI+B(I,K)*W(K)
    1 WX(I)=WXI
      RETURN
      end subroutine
      end module