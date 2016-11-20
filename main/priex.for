      module MOD_PRIEX
      contains
**********  MODULNAME: PRIEX     ******* 24/03/87  19.46.17.******    24 KARTEN
      SUBROUTINE PRIEX (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP)

      IMPLICIT REAL*8(A-H,O-Z) 

      DIMENSION RNE(ND),POPNUM(ND,N)
      CHARACTER MODHEAD*104,LEVEL(N)*10
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A,  20X,'JOB NO.',I5,
     $ //,20X,'RELATIVE NON-LTE POPULATION NUMBERS EXTRAPOLATED ',
     $ 'FROM THE LAST THREE ITERATIONS',/,20X,79('-'))
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2, (LEVEL(J),J=J1,J2)
    2 FORMAT (//,'  L EL.DENS.',10(2X,A10))
      IF (N.LT.J2) PRINT 11
   11 FORMAT (1X)
      PRINT 11
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      PRINT 9,L,RNE(L),(LOG10(POPNUM(L,J)),J=J1,J2)
    9 FORMAT (I3,F7.3,2X,10F12.2)
    3 CONTINUE
      IF (J2.EQ.N) RETURN
      J1=J1+10
      GOTO 4
      END subroutine
      end module