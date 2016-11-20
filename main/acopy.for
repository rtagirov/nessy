      module MOD_ACOPY
      contains
C**********  MODULNAME: ACOPY     ******* 3/07/92  22.14.04.******    11 KARTEN
      SUBROUTINE ACOPY (A1,A2,N,NDIM)
C***  COPY OF ARRAY A1=A2

      implicit real*8(a-h,o-z)

      DIMENSION A1(NDIM,NDIM),A2(NDIM,NDIM)

      DO 1 J=1,NDIM
      DO 2 I=1,NDIM
      A1(J,I)=A2(J,I)
    2 CONTINUE
    1 CONTINUE

      RETURN
      END subroutine
      end module