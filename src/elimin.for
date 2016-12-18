      module MOD_ELIMIN

      contains

      SUBROUTINE ELIMIN(XLAM,EMFLUX,FLUXIN,U,Z,A,B,C,W,BX,WX,XJC,RADIUS,P,BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

!     CALLED BY FORMAL, ETL, FIOSS8, WRCONT
!     FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY
!     TAPE7 = MASS STORAGE FILE FOR FEAUTRIER MATRICES
!     LAST PARAMETER -1 IN WRITMS IS VERY IMPORTANT
!     OTHERWISE "MASS STORAGE LIMIT" EXCEEDED BECAUSE OLD RECORDS ARE NOT OVERWRITTEN.
!     ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY

      use MOD_SETUP
      use MOMENTS
      use MATOPER

      IMPLICIT NONE

      !global variables, intent(inout|out)
      integer ::                  ND, NP
      real*8  ::                  BX(NP, NP, ND)
      real*8  ::                  FLUXIN, WX(NP, ND), EMFLUX
      real*8,dimension(ND, NP) :: U
      real*8,dimension(3, ND)  :: EDDI

      !global variables, intent(in)
      real*8 ::               B(NP, NP), BCORE, DBDR
      real*8,dimension(ND) :: ETA, RADIUS, XJC, OPA, THOMSON
      real*8 ::               P(NP), XLAM, Z(ND,NP)
      real*8,dimension(*) ::  W, C
      real*8,dimension(:) ::  A

      !local variables
      real*8 ::           CORFAC, FL, FLP, H, HPLUS, RL, RLP, RRQ, XK
      integer::           J, JC, JMAX, L, NC2, i

      real*8,parameter :: ONE = 1.D+0, TWO = 2.D0, THREE = 3.D0
     
      do L = 1, ND

         print*, 'elimin XJC here:', L, XJC(L)

      enddo

!      stop

C***  GAUSS-ELIMINATION
      DO L = 1, ND

        CALL SETUP(L,A,B,C,W,JMAX,ND,NP,OPA,ETA,THOMSON,Z,RADIUS,BCORE,DBDR)

        print*, 'elimin WX, check 0:', WX(1, L)

        if (L .NE. 1) then

            CALL MDMV(A, BX(1, 1, L), JMAX, NP)
            CALL MSUB(B, BX(1, 1, L), JMAX, NP)

            print*, 'elimin WX, check 1:', WX(1, L)

            CALL MDV (A, WX(1, L),    JMAX)
            CALL VADD(W, WX(1, L),    JMAX)

            print*, 'elimin WX, check 2:', WX(1, L)

        endif

        CALL INV(JMAX, B(1 : jmax, 1 : jmax))
        CALL MVV (WX(1,L),B,W,JMAX,JMAX,NP)

        print*, 'elimin XJC here 2:', l, XJC(l)

        IF (L.EQ.ND) GOTO 2
        CALL MVMD (BX(1,1,L),B,C,JMAX,JMAX-1,NP)
        !***  COMPRESSING THE MATRIX BX  AND VECTOR WX  INTO THE RANGE OF B AND C
        DO 7 J=1,JMAX
        WX(J,L+1)=WX(J,L)
        DO 7 JC=1,JMAX
  7     BX(JC,J,L+1)=BX(JC,J,L)

      ENDDO
     
C***  BACK SUBSTITUTION
C***  RECENT WX IS THE FEAUTRIER-INTENSITY U AT THE INNER BOUNDARY
    2 CALL MOMENT0(ND,RADIUS,ND,JMAX,Z,WX(1,ND),XJC(ND),.FALSE.)

      do L = 1, ND

         print*, 'elimin XJC here 3:', l, XJC(l)

      enddo

      CALL MOMENT1(RADIUS(ND),JMAX,P,WX(1,ND),H)
      HPLUS=BCORE/two + DBDR/three/OPA(ND)
      CALL MOMENT2(RADIUS(ND), JMAX, P(1 : JMAX), WX(1, ND), XK)

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







      CALL MVV (W,BX(1,1,L),A,JMAX,JMAX-1,NP)
      CALL VADD (WX(1,L),W,JMAX)
C***  WX(J) IS THE FEAUTRIER-INTENSITY U AT RADIUS R(L)
      U(L,:JMAX)=WX(:JMAX,L)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,WX(1,L),XJC(L),.FALSE.)
      CALL MOMENT2(RL, JMAX, P(1 : JMAX), WX(1, L), XK)
!      CALL MOMENT2(RL, JMAX, P(1 : JMAX), U(L, 1 : JMAX), XK)
      EDDI(1, L) = XK / XJC(L)

!      if (isnan(eddi(1, L))) then

      write(*, '(I4,3(2x,e15.7))') l, eddi(1, l), xk, xjc(l)

!         stop 'elimin eddi(1, l) is nan'

!      endif

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
