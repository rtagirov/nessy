      module MOD_PRIPOP

      contains

      SUBROUTINE PRIPOP(LSPOP,WEIGHT,NCHARG,NOM,
     $                  ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)

C***  OUTPUT OF RELATIVE POPULATION NUMBERS AND DEPARTURE COEFFICIENTS
c***  CALLED BY STEAL
     
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in) :: JOBNUM
      DIMENSION POPNUM(ND,N),DEPART(ND,N)
      DIMENSION RNE(ND),ITNE(ND)
      DIMENSION WEIGHT(N),NCHARG(N),NOM(N)
      CHARACTER LEVEL(N)*10
      CHARACTER MODHEAD*104
      CHARACTER PRILINE*130,NUMBER*16
     
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A104  ,20X,'JOB NO.',I5,
     $           //,20X,'RELATIVE NON-LTE POPULATION NUMBERS (LOG) AND',
     $ ' DEPARTURE COEFFICIENTS',/,20X,68('-'))
     
      J1=1
    4 J2=MIN0(N,J1+4)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,4X,'EL.   IT',4X,5(7X,A10,6X))
      IF(N.LT.J2) PRINT 11
   11 FORMAT(1X)
      PRINT 11
     
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      write (PRILINE,9) L,RNE(L),ITNE(L)
    9 FORMAT(I3,F7.3,I3)
     
      DO 5 J=J1,J2
      IF (POPNUM(L,J) .NE. .0) THEN
            write (number,8) LOG10(ABS(POPNUM(L,J))),DEPART(L,J)
    8       FORMAT (F6.2,G10.3)
            ELSE
            NUMBER='****  ZERO  ****'
            ENDIF
      I=18+(J-J1)*23
      PRILINE(I+4:I+19)=NUMBER
      IF (POPNUM(L,J) .LE. .0) THEN
            PRILINE(I:I+3)='NEG'
            ENDIF
C***  LASER WARNING BETWEEN BOUND LEVELS
      IF (J .GT. 1) THEN
      IF ((NCHARG(J-1) .EQ. NCHARG(J)) .AND. (NOM(J-1) .EQ. NOM(J))
     $    .AND. (WEIGHT(J)*POPNUM(L,J-1) .LT. WEIGHT(J-1)*POPNUM(L,J)))
     $      PRILINE(I-2:I-2)='*'
      ENDIF
    5 CONTINUE
     
      PRINT 6,PRILINE
    6 FORMAT (A)
    3 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
     
      IF(J2.EQ.N) GOTO 7
      J1=J1+5
      GOTO 4
    7 CONTINUE
     
      RETURN

      END subroutine

      end module
