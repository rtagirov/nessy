      module MOD_PGRID
      contains
C**********  MODULNAME: PGRID     ******* 24/03/87  21.22.04.******    17 KARTEN
      SUBROUTINE PGRID (NPDIM,NP,ND,R,P)
C***  GRID OF IMPACT-PARAMETER POINTS

      IMPLICIT REAL*8(A-H,O-Z)
      real*8,intent(inout):: P
      integer,intent(out):: NP
      integer,intent(in):: NPDIM,ND
      real*8,intent(in):: R
      DIMENSION  P(NPDIM),R(ND)
C***  NC = NUMBER OF CORE-INTERSECTING RAYS
c      NC=4 (old-value)
      NC=13
!      NC=14
      NP=ND+NC
C***************************************** 
      IF (NP.GT.NPDIM) STOP "PGRID: TOO MANY IMPACT-PARAMETER POINTS"
C***  CORE RAYS EQUELLY SPACED
      D=.9999/FLOAT(NC-1)
      P(1)=0.d0
      P(2)=0.1
      P(3)=0.2
      P(4)=0.3
      P(5)=0.4
      P(6)=0.5
      P(7)=0.6
      P(8)=0.7
      P(9)=0.8
      P(10)=0.91652
      P(11)=0.97980
      P(12)=0.99499
      P(13)=0.99875
!      P(14) = 1.0D0

c      DO 1 J=1,NC
c    1    P(J)=(J-1)*D
c       PRINT *,'J=',J,'!!!!P(J)=',P(J)
      DO 2 L=1,ND
        J=NP+1-L
    2   P(J)=R(L)
c     PRINT *,'L=',L,'J=',J,'!!!!P(J)=',P(J)
      RETURN
      END SUBROUTINE
      end module
