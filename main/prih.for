      module MOD_PRIH
      contains
      SUBROUTINE PRIH (LPRIH,ND,RADIUS,HTOT,TEFF,
     $                  TNEW,TAUROSS,JOBNUM,MODHEAD)
     
      IMPLICIT REAL*8(A-H,O-Z)
 
      DIMENSION RADIUS(ND),HTOT(ND),TNEW(ND),TAUROSS(ND)
      CHARACTER MODHEAD*104
      DATA PI,PISIG/3.141592654d0,5.5411d+4/

      IF (TEFF.LE.0..AND.HTOT(1).GT.0.)
     &             TEFF=(4.*HTOT(1)*RADIUS(1)*RADIUS(1)*PISIG)**0.25
      PRINT 10,MODHEAD,JOBNUM, TEFF
 10   FORMAT (10H1$$$$$$$$$,/,1X,  A,  20X,'JOB NO.',I5,
     $     ///,20X,'ASTROPHYSICAL FLUX AS A FUNCTION OF DEPTH',/,
     $         20X,'=========================================',
     $      //,20X,'TEFF = ',F20.0
     $   //,5X,'DEPTH INDEX  TAU-ROSSELAND  EL. TEMPERATURE  FLUX ',
     $         '(ERG/CM2/S)    T (KELVIN)            T/TEFF'/)

      DO 1 L=1,ND
      IF(((L-1)/LPRIH)*LPRIH.NE.(L-1) .AND. L.NE.ND) GOTO 1
      TL=(4.*HTOT(L)*RADIUS(L)*RADIUS(L)*PISIG)**0.25
      PRINT 13,L,TAUROSS(L),TNEW(L),4.*HTOT(L),TL,TL/TEFF
   13 FORMAT (5X,I10,1PE14.3,0PF14.0,1PE20.5,0PF15.0,F20.5)
    1 CONTINUE
      RETURN
      END subroutine
      end module
