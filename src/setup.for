      module MOD_SETUP

      contains

      SUBROUTINE SETUP(L,A,B,C,W,JMAX,ND,NP,OPA,ETA,THOMSON,Z,RADIUS,BCORE,DBDR)

C***  FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY:
C***  SET UP THE MATRICES  A (DIAGONAL), B (FULL), C (DIAGONAL) AND W (VECTOR)

      use MOD_MOMENT0

      implicit none

      integer,intent(in) :: L, ND, NP

      !inout

      integer,intent(inout) ::                    JMAX
      real*8, intent(inout), dimension(NP) ::     A, C, W
      real*8, intent(inout), dimension(NP, NP) :: B

      !in
      real*8, intent(in) ::                    BCORE, DBDR
      real*8, intent(in), dimension(ND) ::     ETA, OPA, RADIUS, THOMSON
      real*8, intent(in), dimension(ND, NP) :: Z

      !local
      integer :: J, JMM, JS
      real*8 ::  CORFAC, DT, DTP, DTM, DUMMY, ETAL, G
      real*8 ::  PLUSI, WJG, X, XM, XP, ZLJ, ZLMIN, ZLPLUS

      real*8, parameter :: one = 1.D+0, two = 2.d0, three = 3.D0
     
      JMAX = NP + 1 - L
      JMM =  JMAX - 1
     
C***  EVERY L = 1 ... ND
      X=OPA(L)
      G=-X*THOMSON(L)
      ETAL=ETA(L)
     
C***  MEAN INTENSITY INTEGRATION WEIGHTS FROM SUBROUTINE MOMENT0 (VEKTOR W)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,W,DUMMY,.TRUE.)
      DO 1 J=1,JMAX
      WJG=W(J)*G
      DO 1 JS=1,JMAX
    1 B(JS,J)=WJG
      DO 3 J=1,JMAX
    3 W(J)=ETAL
     
      IF(L.EQ.1) GOTO 9
      IF(L.EQ.ND) GOTO 10
     
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1
      XP=(X+OPA(L+1))/two
      XM=(X+OPA(L-1))/two
      DO 2 J=1,JMM
      ZLPLUS=Z(L+1,J)
      ZLJ=Z(L,J)
      ZLMIN=Z(L-1,J)
      DT=two/(ZLMIN-ZLPLUS)
      DTM=XM*(ZLMIN-ZLJ)
      DTP=XP*(ZLJ-ZLPLUS)
      A(J)=DT/DTM
      C(J)=DT/DTP
    2 B(J,J)=B(J,J)+A(J)+C(J)+X
     
C     LAST ROW OF BLOCK, J=JMAX
      ZLMIN=Z(L-1,JMAX)
      DT=ZLMIN*XM
      A(JMAX)=two*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      RETURN
     
C***  OUTER BOUNDARY CONDITION     L = 1
    9 XP=(X+OPA(2))/two	
C***  NONZERO INCIDENT RADIATION
CCCCCC      CORFAC=one-EXP(-RADIUS(1)*X)
C***  AUSSSER BETRIEB ----------------------------- !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CORFAC=.0
C***  WENN C GELOCHT, WIEDER IN BETRIEB !!!  ---------------------------
      DO 8 J=1,JMM
      ZLPLUS=Z(2,J)
      ZLJ=Z(1,J)
      DT=XP*(ZLJ-ZLPLUS)
C***  MODIFICATION FOR NONZERO INCIDENT RADIATION  FROM TRUNCATED LAYERS
      W(J)=ETAL*(one+two*CORFAC/DT)
      C(J)=two*X/DT/DT
      B(JMAX,J)=.0
    8 B(J,J)=B(J,J)+C(J)+X+two*X/DT
      B(JMAX,JMAX)=X
      W(JMAX)=ETAL*CORFAC
      RETURN
     
C***  INNER BOUNDARY CONDITION    L = ND
   10 XM=(X+OPA(ND-1))/two
      DO 14 J=1,JMM
      ZLMIN=Z(ND-1,J)
      ZLJ=Z(ND,J)
      DT=XM*(ZLMIN-ZLJ)
      A(J)=two*X/DT/DT
      B(J,J)=B(J,J)+A(J)+X+two*X/DT
      PLUSI=BCORE+DBDR*ZLJ/X
   14 W(J)=ETAL+PLUSI*two*X/DT
      A(JMAX)=two*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      W(JMAX)=ETAL
      RETURN
      end subroutine
      end module
