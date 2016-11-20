      module MOD_VMALV
      contains
C**********  MODULNAME: VMALV     ******* 03/06/87  12.12.51.******    12 KARTEN
      SUBROUTINE VMALV (VA,VB,V,Q,LMAX)
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
      implicit real*8(a-h,o-z)

      DIMENSION VA(LMAX),VB(LMAX),V(LMAX),Q(LMAX)
      LZ=LMAX-1
      Q(1)=VB(1)*V(1)
      DO 1 L=2,LZ
      Q(L)=VA(L)*V(L-1)+VB(L)*V(L)
    1 CONTINUE
      Q(LMAX)=VA(LMAX)*V(LZ)
      RETURN
      END subroutine
      end module
