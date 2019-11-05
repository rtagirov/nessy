      module MOD_PRIINT

      contains

      SUBROUTINE PRIINT(xjc,xjc2,EDDI,EDDARR,R,ND,XLAMBDA,NF,LSTEP,JOBNUM,MODHEAD)

      use MOD_TRADFUN

      use storextr

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION xjc(ND),XLAMBDA(NF),EDDI(3,ND),R(ND)
      DIMENSION xjc2(ND,NF),EDDARR(3,ND,NF)
      CHARACTER MODHEAD*104

      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,
     $            //,10X,'MEAN INTENSITY',/,10X,14('-'),/,
     $ ' FREQUENCY     DEPTH       J-NUE          T-RAD   ',
     $ '         EDDINGTON     SPHERICITY      EDDINGTON',/,
     $ '   INDEX       INDEX     (ERG/CM+2)       (KELVIN)',
     $ '         FACTOR F      FACTOR  R.R.Q    FACTOR H/J',/)

      DO 4 K=1,NF
c      ENCODE (8,FMT_KEY,NAME) 'XJC ',K
cc    3 FORMAT (A4,I4)
c      CALL READMS (3,XJC,ND,NAME)
c      ENCODE (7,3,NAME) 'EDDI',K
c      CALL READMS (3,EDDI,3*ND,NAME)
      call extrxjc (xjc2,xjc,EDDARR,EDDI,nd,nf,K)
      DO 1 L=1,ND
      IF (((L-1)/LSTEP)*LSTEP .NE. (L-1) .AND. L.NE.ND) GOTO 1
      TRAD=TRADFUN (XLAMBDA(K),XJC(L))
      IF (L.EQ.ND .OR. L.EQ.1) THEN
      PRINT 5,K,L,xjc(L),TRAD,EDDI(1,L),EDDI(2,L)*R(L)*R(L),EDDI(3,L)
    5 FORMAT (2I10,E15.3,F15.0,F15.3,F15.3,F15.3)
      ELSE
      PRINT 5,K,L,xjc(L),TRAD,EDDI(1,L),EDDI(2,L)*R(L)*R(L)
      ENDIF
    1 CONTINUE
    4 PRINT 6
    6 FORMAT (1X)

      return

      END subroutine

      end module
