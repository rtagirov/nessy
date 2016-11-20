       MODULE MOD_INTPL
       contains
       SUBROUTINE INTPL(Y,X,X1,X2,Y1,Y2)
C**************************************
C     CALLED BY INPLCS, CSTABREAD
C     WRITTEN BY MARGIT HABERREITER
C     FOR LINEAR INTERPOLATION
C     INPUT: X1,X2,Y1,Y2,X
C     OUTPUT: Y
C     SLOPE: SLP
C**************************************
!	IMPLICIT REAL*8(a-h,o-z)
	IMPLICIT NONE
	REAL*8,intent(in)::X1,X2,Y1,Y2,X
	REAL*8,intent(out)::Y
	real*8::SLP
	Y=0.
	SLP=(Y2-Y1)/(X2-X1)
	Y=Y1+(SLP*(X-X1))	

	END SUBROUTINE
        END MODULE
