      module MOD_REDCOR
      contains
C**********  MODULNAME: REDCOR    ******* 06/08/87  20.15.06.******    19 KARTEN
      SUBROUTINE REDCOR (POPNUM,POP1,ND,N,RNE,NCHARG,REDUCE)
C*******************************************************************************
C***  REDUCTION OF THE CORRECTIONS OF POP-NUMBERS BY FACTOR "REDUCE"
C******************************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION POPNUM(ND,N),POP1(ND,N),RNE(ND),NCHARG(N)
     
      DO 2 L=1,ND
      RNEL=.0
      DO 1 J=1,N
      POPNUM(L,J)=REDUCE*POPNUM(L,J)+(1.-REDUCE)*POP1(L,J)
      RNEL=RNEL+NCHARG(J)*POPNUM(L,J)
    1 CONTINUE
      RNE(L)=RNEL
    2 CONTINUE
     
      RETURN
      END subroutine
      end module