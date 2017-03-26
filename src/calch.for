      module MOD_CALCH

      contains

      SUBROUTINE CALCH(ND,NP,OPA,Z,P,U,VL,VJL,RADIUS,HNU,FLUXIN,EMFLUX)

      use MOMENTS

      IMPLICIT REAL*8(A-H,O-Z)

C***
      parameter (two=2.d0, four=4.d0)
C***
      DIMENSION RADIUS(ND),OPA(ND)
      DIMENSION Z(ND,NP),P(NP),U(ND,NP),HNU(ND),VL(NP),VJL(NP,ND)

C***  EVERY L = 1 ... ND
      DO 1 L=1,ND
      JMAX=NP+1-L
      JMM=JMAX-1

      IF(L.EQ.1) THEN
      HNU(1)=EMFLUX/four/RADIUS(1)/RADIUS(1)
      ELSE IF (L.GT.1 .AND. L.LT.ND) THEN
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1
      XP=(OPA(L)+OPA(L+1))/two
      XM=(OPA(L)+OPA(L-1))/two
      DO 2 J=1,JMM
      ZLPLUS=Z(L+1,J)
      ZLJ=Z(L,J)
      ZLMIN=Z(L-1,J)
      DTM=XM*(ZLMIN-ZLJ)
      DTP=XP*(ZLPLUS-ZLJ)
      DUM=U(L-1,J)-U(L,J)
      DUP=U(L+1,J)-U(L,J)
      DDT=DTP-DTM
      DDU=DUM*DTP/DTM-DUP*DTM/DTP
c*** Mihalas p152 6-15 du/dtau=v
      VL(J)=-DDU/DDT
    2 CONTINUE
      VL(JMAX)=0.d0
      DO 3 J=1,JMAX
    3 VJL(J,L)=VL(J)  
C*** CALCULATE H = INTEG 0 TO 1  V P/R/R DP
      CALL MOMENT1 (RADIUS(L),JMAX,P,VL,HNU(L))

      ELSE IF(L.EQ.ND) THEN
      HNU(ND)=FLUXIN/four
      ENDIF
    1 CONTINUE
C
      RETURN
      end subroutine
      end module
