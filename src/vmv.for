      module MOD_VMV
      contains
C**********  MODULNAME: VMV     ******* 24/06/92  22.14.04.******    11 KARTEN
      SUBROUTINE VMV (VSUM,V1,V2,N,NDIM)
C***  MULTIPLICATION RESULT = VECTOR * VEKTOR  --  SUM = V1 * V2^T
      implicit real*8(a-h,o-z)

      DIMENSION V1(NDIM),V2(NDIM)
      VSUM=.0
      DO 1 I=1,N
      VSUM=VSUM+V1(I)*V2(I)
    1 CONTINUE
      RETURN
      END subroutine
      end module