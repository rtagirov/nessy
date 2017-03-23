      module MOD_PRIINTL

      contains

      SUBROUTINE PRIINTL(N,LEVEL,WEIGHT,EINST,LASTIND,LINE,NLINE,
     $                   INDLOW,INDNUP,ELEVEL,ND,XJL,LSINT,JOBNUM,MODHEAD)
C***********************************************************************
C***  PRINTOUT OF THE LINE INTENSITIES
C***********************************************************************
      use MOD_TRADFUN

      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION WEIGHT(N),EINST(N,N),ELEVEL(N)
      DIMENSION XJL(ND,lastind)
c      ,LINE(1)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      CHARACTER MODHEAD*104,LEVEL(N)*10
      character*7 LINE(LASTIND), name
     
      PRINT 4,MODHEAD,JOBNUM
    4 FORMAT (10H1$$$$$$$$$,/,1X,A  ,20X,'JOB NO.',I5,
     $ //,30X,'MEAN INTENSITIES IN THE LINES',/,30X,29('='))
     
C***  LOOP FOR EVERY LINE TRANSITION  -----------------------------------
      DO 16 IND=1,LASTIND
    1 FORMAT (A4,I3)
c      ENCODE (7,1,NAME) 'LINE',IND
     	write (name,1) 'LINE',ind
     
C***  LINE IS PRINTED ONLY IF IT IS QUOTED ON AN INPUT OPTION CARD
      DO 7 NLHELP=1,NLINE
      NL=NLHELP
      IF (LINE(NLHELP) .EQ. NAME) GOTO 6
    7 CONTINUE
      GOTO 16
    6 CONTINUE
     
      I=INDLOW(IND)
      J=INDNUP(IND)
C***  NO PRINTOUT FOR RUDIMENTAL LINES
      IF (EINST(I,J) .EQ. -2.) GOTO 16
     
      XLAM=1.E8/(ELEVEL(J)-ELEVEL(I))
      F = 1.499E-16*XLAM*XLAM*EINST(J,I)*WEIGHT(J)/WEIGHT(I)
c      ENCODE (7,1,NAME) 'XJL ',IND
c      CALL READMS (3,XJL,ND,NAME)
      PRINT 5, LINE(NL),XLAM,LEVEL(I),LEVEL(J),EINST(J,I),F
    5 FORMAT (//,10X,A10,5X,'LAMBDA =',F10.2,' (ANGSTROEM)',
     $ 5X,'FROM LEVEL ',A10,' (LOW)  TO LEVEL ',A10,' (UP)',
     $ /,10X,103('-'),/,
     $ 17X,'EINSTEIN-COEFFICIENT   A(UP-LOW) =',1P,E12.4,0P,
     $ ' (PER SECOND)', 7X,'OSCILLATOR STRENGTH  F =',F6.3,//,
     $ 20X,'DEPTH         J-NUE         TRAD',/,
     $ 20X,'INDEX       (ERG/CM+2)      (KELVIN)',/)
     
      DO 10 L=1,ND
      IF (((L-1)/LSINT)*LSINT .NE. (L-1) .AND. L.NE.ND) GOTO 10
      TRAD=TRADFUN (XLAM,XJL(L,ind))
      PRINT 11,L,XJL(L,ind),TRAD
   11 FORMAT (20X,I3,E18.3,F14.0)
   10 CONTINUE
     
   16 CONTINUE
C***  ENDLOOP  ----------------------------------------------------------
     
      RETURN
      END subroutine
      end module
