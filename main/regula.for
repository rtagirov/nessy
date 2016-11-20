      module MOD_REGULA
      contains
C**********  MODULNAME: REGULA    ******* 24/03/87  22.00.45.******    46 KARTEN
      SUBROUTINE REGULA(F,X,Y,X1,X2,EPS)
C***********************************************************************
C***  THIS ROUTINE CALCULATES THE SOLUTION X OF F(X)=Y IN THE INTERVAL
C***  (X1,X2) | METHOD:
C***  -------- REGULA FALSI --------
C***  USING BISECTION STEPS TO GUARANTEE CONVERGENCE, PRECISION IN X: EPS
C***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
 
      DATA ANULL/0.d0/
      LOGICAL BI
     
      BI=.TRUE.
      A=X1
      B=X2
      FA=F(A)-Y
      IF(FA.EQ..0) GOTO 3
      FB=F(B)-Y
      IF(FB.EQ..0) GOTO 4
      IF (FA*FB.GT..0) THEN
            PRINT 10,A,B,FA,FB
 10         FORMAT(//,10X,20(1H*),' ERROR IN REGULA ',20(1H*),/,
     $      10X,'INTERVAL:      A =',E15.5,10X,'   B =',E15.5,/,
     $      10X,'FUNCTION:   F(A) =',E15.5,10X,'F(B) =',E15.5,/,
     $      10X,20(1H*),' ERROR IN REGULA ',20(1H*),/)
            fnull = F(ANULL)
            fa=f(a)
            fb=f(b)
            print *,' Y    = ',Y
            print *,' F(0) = ',fnull
            print *,' F(A) = ',fa
            print *,' F(B) = ',fb

c** creat run time error to find out which was the calling routine
            print *,' run time error deliberate!'
            error=FB/ANULL
            error=error**2.
            print *,error
            STOP 'ERROR'
            ENDIF
    1 D=A-B
      IF(ABS(D).LT.EPS) GOTO 5
      BI=.NOT.BI
      IF(BI) X=A-FA*D/(FA-FB)
      IF(.NOT.BI) X=.5*(A+B)
      FX=F(X)-Y
      IF(FX.EQ..0) RETURN
      IF(FX*FA.GT..0) GOTO 2
      B=X
      FB=FX
      GOTO 1
    2 A=X
      FA=FX
      GOTO 1
    3 X=A
      RETURN
    4 X=B
      RETURN
    5 X=(A+B)*.5
      RETURN
      END subroutine
      end module