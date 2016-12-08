      module MOD_PRICC
      contains

      SUBROUTINE PRICC (ND,NF,WCHARM,DELTAC,MODHEAD,JOBNUM)
C***  PRINTOUT OF SCHARMER CONTINUUM CORES (WEIGHT FUNCTION)

	IMPLICIT NONE
	CHARACTER,intent(in) :: MODHEAD*104
	INTEGER,intent(in)   :: JOBNUM, NF, ND
      INTEGER :: K, L
	REAL*8,intent(in)   :: DELTAC
      REAL*8,dimension(ND,NF),intent(in)    :: WCHARM
      
     
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,//)
      PRINT 2,DELTAC, NF
    2 FORMAT (10X,'CONTINUUM SCHARMER WEIGHTS (FIRST DECIMAL DIGIT) -',
     $      '   DELTAC=',F5.1,' NF=',I5//)
      PRINT 21,(K/10,K=1,NF)
      PRINT 21,(K-K/10*10,K=1,NF)
   21 FORMAT (10X,100I1)
      PRINT 23
   23 FORMAT (1X)
      DO L=1,ND
      	PRINT '(I5,5X,100I1)',L,(INT(WCHARM(L,K)*10.),K=1,NF)
      ENDDO
      RETURN
      END subroutine
      end module
