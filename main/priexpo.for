      module MOD_PRIEXPO
      contains
C**********  MODULNAME: PRIEXPO   ******* 24/03/87  21.29.35.******    97 KARTEN
      SUBROUTINE PRIEXPO (POPNUM,POP1,POP2,LEVEL,N,ND,MODHEAD,JOBNUM,
     $                    LSEXPO )
C***  PRINTOUT OF RELATIVE CORRECTIONS WHICH WOULD RESULT IF THE CORRECTIONS
C***  FOLLOW A GEOMETRICAL SERIES

      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in) :: JOBNUM
      DIMENSION POPNUM(ND,N),POP1(ND,N),POP2(ND,N)
      CHARACTER LEVEL(N)*10
      CHARACTER MODHEAD*104
      CHARACTER PRILINE*130,NUMBER*12,LEVMAX*10,NUMMAX*12
     
C***  THIS SUBROUTINE MAY NOT WORK AS LONG AS JOBNUM .LT. 10
      IF (JOBNUM .LT. 10) RETURN
     
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,A,20X,'JOB NO.',I5,//,10X,
     $   'ERROR ESTIMATE: RELATIVE CORRECTIONS AS EXPECTED FROM AITKEN',
     $      ' EXTRAPOLATION',/,10X,80('-'))
     
C*** POP2 IS READ (AND SHIFTED TO POP2) AT THE START OF steal
C      CALL READMS (3,POP2,ND*N,4HPOP2 )
     
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,' DEPTH',10(2X,A10))
      IF (N .LT. J2) PRINT 5
    5 FORMAT (1X)
      PRINT 5
     
      CORMAX=.0
      LMAX=0
      DO 3 L=1,ND
      ENCODE (130,6,PRILINE) L
    6 FORMAT (I6)
      DO 12 J=J1,J2
      IF (POPNUM(L,J) .LE. .0) GOTO 23
      D1=POPNUM(L,J)-POP1(L,J)
      D2=POP1  (L,J)-POP2(L,J)
      IF (D2 .EQ. .0) THEN
            IF (D1 .EQ. .0) GOTO 20
            IF (D1 .GT. .0) GOTO 21
            GOTO 22
            ENDIF
      Q=D1/D2
      IF (Q .GE. 1. ) THEN
            IF (D1 .GT. .0) GOTO 21
            GOTO 22
            ENDIF
      COR=D1/POPNUM(L,J)/(1.-Q)
      ENCODE (12,11,NUMBER) COR
   11 FORMAT (F12.4)
      IF (NUMBER .EQ. '      0.0000') NUMBER='      0     '
      IF (ABS(COR) .GT. CORMAX) THEN
            CORMAX=ABS(COR)
            NUMMAX=NUMBER
            LEVMAX=LEVEL(J)
            LMAX=L
            ENDIF
      GOTO 16
     
   20 NUMBER='           0'
      GOTO 16
     
   21 NUMBER='     +I     '
      GOTO 17
     
   22 NUMBER='     -I     '
      GOTO 17
     
   23 NUMBER='   UNDEFINED'
     
   17 CORMAX=1.d99
      LMAX=L
      NUMMAX=NUMBER
      LEVMAX=LEVEL(J)
     
   16 CONTINUE
      I=7+(J-J1)*12
      PRILINE(I:I+11)=NUMBER
   12 CONTINUE
     
      IF (((L-1)/LSEXPO)*LSEXPO .EQ. (L-1) .OR.  L .EQ. ND )
     $PRINT 13,PRILINE
   13 FORMAT (A)
    3 CONTINUE
     
      IF (J2 .LT. N) THEN
            J1=J1+10
            GOTO 4
            ENDIF
     
C***  PRINTOUT OF MAXIMUM CORRECTION
      PRINT 9,NUMMAX,LEVMAX,LMAX
    9 FORMAT (/,10X,'LARGEST CORRECTION ENCOUNTERED:',A12,
     $   '   FOR LEVEL: ',A10,' AT DEPTH POINT L=',I3)
      RETURN
      END subroutine
      end module
