      module MOD_LPPLOT
      contains
      SUBROUTINE LPPLOT (KANAL,X0,DX,NF,F,KEYWORD,MODHEAD,JOBNUM)
c234567890 234567890 234567890 234567890 234567890 234567890 234567890     
C************************************************************************
C***  LINE PRINTER PLOT OF ONE FUNCTION F(X)
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION F(NF)
      DIMENSION IFORM(8)
c     $,IFORM1(4),IFORM2(4)
      character*28 iform1,iform2
      CHARACTER*1 LP(123)
      CHARACTER MODHEAD*104
	character keyword*80
C***  STATEMENT FUCTION, GIVING THE TABULATOR POSITION OF THE PROFILE POINT
      L1FUN(K)=3+INT(.5+10*K2*F(K)*SCALE)
     
C***  HEADER LINE
      WRITE (KANAL,2) MODHEAD,JOBNUM,KEYWORD
    2 FORMAT (10H1$$$$$$$$$,/,1X,A,20X,'JOB NO.',I5,
     $                      /,50('='),2X,A8 ,2X,50('='),/)
     
C***  RE-SCALING, IF FUNCTION MAXIMUM GT. 3.
      SCALE=1.
      FMAX=1.
      DO 1 K=1,NF
    1 IF (F(K).GT.FMAX) FMAX=F(K)
      IF (FMAX.LT.3.) GOTO 5
      SCALE=1./FMAX
      WRITE (KANAL,3) FMAX
    3 FORMAT (1X,10('-'),' FUNCTION MAXIMUM',1P,E10.3,
     $ ' SCALED TO UNITY ',10('-'),/)
      FMAX=1.
     
C***  FIRST AXIS
    5 K1=1+INT(10*FMAX)
      K2=INT(120./FLOAT(K1))
      K3=K2-1
      N1=10*K2-1
      N2=(K1-10)*K2-1
    4 FORMAT(10H(9X,'0.0',,I3,1HF,I3,3H.1))
      write (iform1,4) K1,K2
c      ENCODE(30,4,IFORM1) K1,K2
    6 FORMAT(9H(10X,'I',,I3,1H(,I3,12H('-'),'I'))))
      write (iform2,6)  K1,K3
c      ENCODE(30,6,IFORM2) K1,K3
      WRITE (KANAL,IFORM1) (I*.1,I=1,K1)
      WRITE (KANAL,IFORM2)
     
      N1=10*K2+3
      N2=K1*K2+3
     
C***  RESTRICT THE PLOT TO THE PART WHERE THE PROFILE DIFFERS FROM UNITY
      DO 20 K=1,NF-1
      KFIRST=K
      L1=L1FUN(K)
      L1NEXT=L1FUN(K+1)
      IF (L1 .NE. N1 .OR. L1NEXT .NE. N1) GOTO 21
   20 CONTINUE
   21 CONTINUE
     
      DO 22 K=NF,KFIRST+1,-1
      KLAST=K
      L1=L1FUN(K)
      L1LAST=L1FUN(K-1)
      IF (L1 .NE. N1 .OR. L1LAST .NE. N1) GOTO 23
   22 CONTINUE
   23 CONTINUE
     
C***  LOOP -------------------------------------------------------------
      DO 8 K=KFIRST,KLAST
      DO 11 L=2,122
   11 LP(L)=' '
      L1=L1FUN(K)
      IF (L1.LT.1) L1=1
      IF (L1.GT.123) L1=123
     
      LP(3)='I'
      LP(N2)='I'
      LP(N1)='.'
      IF (LP(L1).EQ.' ') LP(L1)='*'
      IF (LP(L1).EQ.'I') LP(L1)='X'
      IF (LP(L1).EQ.'.') LP(L1)='I'
    8 WRITE (KANAL,7) X0+K*DX,(LP(I),I=1,123)
    7 FORMAT(F8.1,123A1)
     
C***  SECOND AXIS
      WRITE (KANAL,IFORM2)
      WRITE (KANAL,IFORM1) (I*.1,I=1,K1)
     
      RETURN
      END subroutine
      end module
