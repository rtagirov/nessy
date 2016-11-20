      module MOD_VMT
      contains
C**********  MODULNAME: VMT     ******* 24/06/92  22.14.04.******    11 KARTEN
      SUBROUTINE VMT (V2,A,V1,N,NDIM)
C***  MULTIPLICATION VECTOR = MATRIX (FULL) * VEKTOR TRANSFORMIERT --  V2 = A * V1^T
      implicit real*8(a-h,o-z)

      DIMENSION V1(NDIM),V2(NDIM),A(NDIM,NDIM)
      DO 1 J=1,N
      VSUM=.0
      DO 2 I=1,N
      VSUM=VSUM+V1(I)*A(J,I)
    2 CONTINUE
      V2(J)=VSUM
    1 CONTINUE
      RETURN
      END subroutine
      end module