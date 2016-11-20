      module MOD_PRIMINT
      contains
C**********  MODULNAME: PRIMINT   ******* 24/03/87  21.29.41.******    21 KARTEN
      SUBROUTINE PRIMINT (XJC,ND,XLAMBDA,NF,K,LSINT,EDDI,JOBNUM,MODHEAD)
      use MOD_TRADFUN
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION XJC(ND,NF),XLAMBDA(NF),EDDI(3,ND)
      CHARACTER MODHEAD*104
      PRINT 2,MODHEAD,JOBNUM,XLAMBDA(K)
    2 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,
     $ //,10X,'MEAN INTENSITY FOR LAMBDA =',F15.2/,10X,14('-'),/,
     $ ' FREQUENCY     DEPTH       J-NUE          T-RAD   ',  /,
     $ '   INDEX       INDEX     (ERG/CM+2)       (KELVIN)',/)
C      DO 4 K=1,NF
      DO 1 L=1,ND
      IF (((L-1)/LSINT)*LSINT .NE. (L-1) .AND. L.NE.ND) GOTO 1
      XJCLK=XJC(L,K)
      TRAD=TRADFUN (XLAMBDA(K),XJC L K )
      PRINT 5,K,L,XJCLK ,TRAD , (EDDI(J,L),J=1,3)
    1 CONTINUE
    5 FORMAT (2I10,E15.3,F15.0,3F15.5)
    4 PRINT 6
    6 FORMAT (1X)
      RETURN
      END subroutine
      end module