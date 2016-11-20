      module MOD_VMD
      contains
C**********  MODULNAME: VMD     ******* 24/06/92  22.14.04.******    11 KARTEN
      SUBROUTINE VMD (A,V1,V2,N,NDIM)
C***  MULTIPLICATION MATRIX = VECTOR * VEKTOR, DYADIC PRODUCT -- A = V1^T * V2
      implicit real*8(a-h,o-z)

      DIMENSION V1(NDIM),V2(NDIM),A(NDIM,NDIM)
      DO 1 J=1,N
      DO 2 I=1,N
      A(J,I)=V1(J)*V2(I)
    2 CONTINUE
    1 CONTINUE
      RETURN
      END subroutine
      end module