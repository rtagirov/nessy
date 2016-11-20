      module MOD_CHANGE
      contains
      SUBROUTINE CHANGE (array1,array2,N)
C***  THIS SUBROUTINE copies an array

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION array1(N), array2(N)
      DO 1 I=1,N
         array2(i)=array1(i)
    1 CONTINUE
      RETURN
      END subroutine
      end module