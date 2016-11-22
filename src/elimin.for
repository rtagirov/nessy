      module MOD_ELIMIN
      contains
      SUBROUTINE ELIMIN (XLAM,EMFLUX,FLUXIN,U,Z,
     $          A,B,C,W,BX,WX,XJC,RADIUS,P,BCORE,DBDR,
     $                       OPA,ETA,THOMSON,EDDI,ND,NP,NPDIM)
C***  CALLED BY FORMAL, ETL, FIOSS8, WRCONT
C***  FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY
C***  TAPE7 = MASS STORAGE FILE FOR FEAUTRIER MATRICES
C***  LAST PARAMETER -1 IN WRITMS IS VERY IMPORTANT !
C***  OTHERWISE "MASS STORAGE LIMIT" EXCEEDED BECAUSE OLD RECORDS ARE
C***  NOT OVERWRITTEN .
C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY !
      use MOD_INV
      use MOD_SETUP
      use MOD_MOMENT0
      use MOD_MOMENT1
      use MOD_MOMENT2
      use MOD_MDMV
      use MOD_MSUB
      use MOD_MDV
      use MOD_VADD
      use MOD_MVV
      use MOD_MVMD
      IMPLICIT NONE
      !global variables, intent(inout|out)
      integer :: ND,NP,NPDIM
      real*8  ::BX(NPDIM,NPDIM,ND)
      real*8  ::FLUXIN,WX(NPDIM,ND),EMFLUX
      real*8,dimension(ND,NP) ::U
      real*8,dimension(3,ND)  ::EDDI
      !global variables, intent(in)
      real*8 :: B(:,:),BCORE,DBDR
      real*8,dimension(ND) :: ETA,RADIUS,XJC,OPA,THOMSON
      real*8 :: P(NPDIM),XLAM,Z(ND,NP)
      real*8,dimension(*) :: W,C
      real*8,dimension(:) :: A
      !local variables
      real*8 :: CORFAC,FL,FLP,H,HPLUS,RL, RLP, RRQ, XK
      integer:: J,JC,JMAX,L,NC2
      real*8,parameter :: ONE = 1.D+0, TWO = 2.D0, THREE = 3.D0
     
C***  GAUSS-ELIMINATION
      DO L=1,ND
        CALL SETUP (L,A,B,C,W,JMAX,ND,NP,NPDIM,OPA,ETA,THOMSON,
     $                  XLAM,Z,RADIUS,BCORE,DBDR )
        if (L.NE.1) then
          CALL MDMV (A,BX(1,1,L),JMAX,NPDIM)
          CALL MSUB (B,BX(1,1,L),JMAX,NPDIM)
          CALL MDV (A,WX(1,L),JMAX)
          CALL VADD (W,WX(1,L),JMAX)
        endif
!        CALL INV (JMAX,NPDIM,B)
        CALL INV (JMAX,B)
        CALL MVV (WX(1,L),B,W,JMAX,JMAX,NPDIM)
        IF (L.EQ.ND) GOTO 2
        CALL MVMD (BX(1,1,L),B,C,JMAX,JMAX-1,NPDIM)
        !***  COMPRESSING THE MATRIX BX  AND VECTOR WX  INTO THE RANGE OF B AND C
        DO 7 J=1,JMAX
        WX(J,L+1)=WX(J,L)
        DO 7 JC=1,JMAX
  7     BX(JC,J,L+1)=BX(JC,J,L)
c      DO 7 J=1,JMAX
c      JC=1+(J-1)*JMAX
c    7 CALL EQUAL (B (JC)  ,BX(1,J),JMAX)
c      CALL EQUAL (B (JMAX*JMAX+1)  ,WX,JMAX)
c      LL=JMAX*(JMAX+1)
c      CALL WRITMS (7,B ,LL,L,-1,IDUMMY,IERR)
      ENDDO
     
C***  BACK SUBSTITUTION
C***  RECENT WX IS THE FEAUTRIER-INTENSITY U AT THE INNER BOUNDARY
    2 CALL MOMENT0 (ND,RADIUS,ND,JMAX,Z,WX(1,ND),XJC(ND),.FALSE.)
      CALL MOMENT1 (RADIUS(ND),JMAX,P,WX(1,ND),H)
      HPLUS=BCORE/two + DBDR/three/OPA(ND)
      CALL MOMENT2 (RADIUS(ND),JMAX,P,WX(1,ND),XK)
C***  EDDI(1,L) IS THE EDDINGTON FACTOR  F = K / J
C***  EDDI(2,L) IS THE SPHERICITY FACTOR Q
C***  EDDI(3,L) IS THE EDDINGTON FACTOR H / J  (ONLY AT THE BOUNDARIES)
C***  EDDI(3,ND-1) IS THE OUTWARD FLUX HPLUS AT THE INNER BOUNDARY
      EDDI(1,ND)=XK/XJC(ND)
      EDDI(2,ND)=one
      RRQ=one
      EDDI(3,ND)=H/XJC(ND)
      EDDI(3,ND-1)=HPLUS
      FLUXIN=4*(HPLUS-H)
      FL=three-one/EDDI(1,ND)
      U(ND,:JMAX)=WX(:JMAX,ND)
      A(:JMAX) = WX(:JMAX,ND)
      NC2=NP-ND+2
     
C***  L = ND-1 ... 1
      DO 4 JMAX=NC2,NP
      L=NP+1-JMAX
      RL=RADIUS(L)
c      LL=JMAX*(JMAX+1)
c      CALL READMS (7,B ,LL,L,IERR)
C***  DECOMPRESSING THE MATRIX BX AND VECTOR WX  OUT OF B (AND C)
c      DO 8  J=1,JMAX
c      JC=1+(J-1)*JMAX
c    8 CALL EQUAL (BX(1,J),B(JC)  ,JMAX)
c      CALL EQUAL (WX,B (JMAX*JMAX+1)  ,JMAX)
      CALL MVV (W,BX(1,1,L),A,JMAX,JMAX-1,NPDIM)
      CALL VADD (WX(1,L),W,JMAX)
C***  WX(J) IS THE FEAUTRIER-INTENSITY U AT RADIUS R(L)
      U(L,:JMAX)=WX(:JMAX,L)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,WX(1,L),XJC(L),.FALSE.)
      CALL MOMENT2 (RL,JMAX,P,WX(1,L),XK)
      EDDI(1,L)=XK/XJC(L)
C***  THIS IS AN INGENIOUS (;) RECURSION FORMULA FOR THE SPHERICITY FACTOR !
      RLP=RADIUS(L+1)
      FLP=FL
      FL=three-one/EDDI(1,L)
      RRQ=RRQ *EXP(FL-FLP)*(RL/RLP)**((FLP*RL-FL*RLP)/(RL-RLP))
      EDDI(2,L)=RRQ/RL/RL
      A(:JMAX) = WX(:JMAX,L)
    4 enddo
     
      CALL MOMENT1 (RADIUS(1),NP,P,WX(1,1),H)
C***  MODIFICATION FOR NONZERO INCIDENT RADIATION  FROM TRUNCATED LAYERS
CCCCCC      CORFAC=one-EXP(-RADIUS(1)*OPA(1))
C***  AUSSSER BETRIEB ----------------------------- !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CORFAC=.0
C***  WENN C GELOCHT, WIEDER IN BETRIEB !!!  ---------------------------
      H=H-0.5*CORFAC*ETA(1)/OPA(1)
C***  EMFLUX IS THE EMERGENT FLUX, BUT RELATED TO THE INNER RADIUS
      EMFLUX=4*H*RADIUS(1)*RADIUS(1)
      EDDI(3,1)=H/XJC(1)
     
      RETURN
      END subroutine
      end module
