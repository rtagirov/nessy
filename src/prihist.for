      module MOD_PRIHIST
      contains
C**********  MODULNAME: PRIHIST   ******* 24/03/87  21.29.39.******    24 KARTEN
      SUBROUTINE PRIHIST (MODHEAD,JOBNUM)
C***  PRINTOUT OF THE MODEL HISTORY

      IMPLICIT REAL*8(A-H,O-Z)

      CHARACTER MODHEAD*104, card*120

      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,1X,  A  ,10X,'JOB NO.',I5,//,11X,
     $ 'M O D E L   H I S T O R Y',/,11X,25('='),//)
     
      open (7,file='MODHIST',status='old')
  8   read (7,'(A120)',end=11) card
      print 4, card
	goto 8
 11   continue
      close (7)
    4 FORMAT (11X,A120 )

      RETURN
      END subroutine
      end module
