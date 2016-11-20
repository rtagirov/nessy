      module MOD_VDIFF
      contains
C**********  MODULNAME: VDIFF   ******* 24/06/92  22.14.04.******    11 KARTEN
      SUBROUTINE VDIFF (A1,VDI,N,NDIM)
C***  DIVISION MATRIX = MATRIX / FLOAT  --  A1 = A2 / VDIFF
      implicit real*8(a-h,o-z)

      DIMENSION A1(NDIM,NDIM)
      IF (VDI .EQ. .0) GOTO 3
      DO 2 J=1,N
      DO 1 I=1,N
      A1(J,I)=A1(J,I)/VDI
    1 CONTINUE
    2 CONTINUE
    3 CONTINUE
      RETURN
      END subroutine
      end module