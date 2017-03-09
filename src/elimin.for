      module MOD_ELIMIN

      contains

      SUBROUTINE ELIMIN(XLAM,EMFLUX,FLUXIN,U,Z,XJC,RADIUS,P,BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

!     CALLED BY FORMAL, ETL, FIOSS, WRCONT
!     FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY

      use MOMENTS
      use MATOPER

      implicit none

      integer, intent(in) :: ND, NP

      real*8, intent(in) ::                BCORE, DBDR
      real*8, intent(in) ::                P(NP), XLAM, Z(ND, NP)
      real*8, dimension(ND), intent(in) :: ETA, RADIUS, OPA, THOMSON


      real*8, intent(out) ::                    FLUXIN, EMFLUX

      real*8, dimension(3, ND),  intent(out) :: EDDI

      real*8, dimension(ND, NP), intent(out) :: U

      real*8, dimension(ND),     intent(out) :: XJC

      real*8, allocatable, dimension(:) ::       A, W, C
      real*8, allocatable, dimension(:, :) ::    B, WX
      real*8, allocatable, dimension(:, :, :) :: BX

      real*8  :: FL, FLP, H, HPLUS, RL, RLP, RRQ, XK
      integer :: J, JC, JMAX, L, NC2

      real*8, parameter :: ONE = 1.D+0, TWO = 2.D0, THREE = 3.D0

      allocate(A(NP), C(NP), W(NP))
      allocate(B(NP, NP))

      allocate(BX(NP, NP, ND)); allocate(WX(NP, ND))

      A(1 : NP) = 0.0D0
      C(1 : NP) = 0.0D0
      W(1 : NP) = 0.0D0

      B(1 : NP,  1 : NP) = 0.0D0

      WX(1 : NP, 1 : ND) = 0.0D0

      BX(1 : NP, 1 : NP, 1 : ND) = 0.0D0

!     GAUSS-ELIMINATION
      DO L = 1, ND

         CALL SETUP(L,A,B,C,W,JMAX,ND,NP,OPA,ETA,THOMSON,Z,RADIUS,BCORE,DBDR)

         if (L .NE. 1) then

             CALL MDMV(A(1 : JMAX),           BX(1 : JMAX, 1 : JMAX, L), JMAX)
             CALL MSUB(B(1 : JMAX, 1 : JMAX), BX(1 : JMAX, 1 : JMAX, L), JMAX)

             CALL  MDV(A(1 : JMAX), WX(1 : JMAX, L), JMAX)
             CALL VADD(W(1 : JMAX), WX(1 : JMAX, L), JMAX)

         endif

         CALL INV(JMAX, B(1 : JMAX, 1 : JMAX))

         CALL MVV(WX(1 : JMAX, L), B(1 : JMAX, 1 : JMAX), W(1 : JMAX), JMAX, JMAX)

         IF (L .EQ. ND) GOTO 2

         CALL MVMD(BX(1 : JMAX, 1 : JMAX, L), B(1 : JMAX, 1 : JMAX), C(1 : JMAX), JMAX, JMAX - 1)

!        COMPRESSING THE MATRIX BX  AND VECTOR WX  INTO THE RANGE OF B AND C
         DO J = 1, JMAX

            WX(J, L + 1) = WX(J, L)

           DO JC = 1, JMAX

              BX(JC, J, L + 1) = BX(JC, J, L)

           ENDDO

         ENDDO

      ENDDO

      stop

!     BACK SUBSTITUTION
!     RECENT WX IS THE FEAUTRIER-INTENSITY U AT THE INNER BOUNDARY
    2 CALL MOMENT0_ELIMIN(ND, RADIUS, ND, JMAX, Z(1 : ND, 1 : JMAX), WX(1 : JMAX, ND), XJC(ND))

      CALL MOMENT1(RADIUS(ND), JMAX, P(1 : JMAX), WX(1 : JMAX, ND), H)

      HPLUS = BCORE / two + DBDR / three / OPA(ND)

      CALL MOMENT2(RADIUS(ND), JMAX, P(1 : JMAX), WX(1 : JMAX, ND), XK)

!     EDDI(1,L) IS THE EDDINGTON FACTOR  F = K / J
!     EDDI(2,L) IS THE SPHERICITY FACTOR Q
!     EDDI(3,L) IS THE EDDINGTON FACTOR H / J  (ONLY AT THE BOUNDARIES)
!     EDDI(3,ND-1) IS THE OUTWARD FLUX HPLUS AT THE INNER BOUNDARY

      EDDI(1, ND) = XK / XJC(ND)
      EDDI(2, ND) = one
      EDDI(3, ND) = H / XJC(ND)
      EDDI(3, ND - 1) = HPLUS

      FLUXIN = 4 * (HPLUS - H)
      FL = three - one / EDDI(1, ND)

      U(ND, 1 : JMAX) = WX(1 : JMAX, ND)
      A(1 : JMAX) = WX(1 : JMAX, ND)

      NC2 = NP - ND + 2
     
      RRQ = one

!     L = ND-1 ... 1
      DO JMAX = NC2, NP

         L = NP + 1 - JMAX

         RL = RADIUS(L)

         CALL MVV(W(1 : JMAX), BX(1 : JMAX, 1 : JMAX, L), A(1 : JMAX), JMAX, JMAX - 1)
         CALL VADD(WX(1 : JMAX, L), W(1 : JMAX), JMAX)

!        WX(J) IS THE FEAUTRIER-INTENSITY U AT RADIUS R(L)
         U(L, 1 : JMAX) = WX(1 : JMAX, L)

         CALL MOMENT0_ELIMIN(ND, RADIUS, L, JMAX, Z(1 : ND, 1 : JMAX), WX(1 : JMAX, L), XJC(L))

         CALL MOMENT2(RL, JMAX, P(1 : JMAX), WX(1 : JMAX, L), XK)

         EDDI(1, L) = XK / XJC(L)

         if (isnan(eddi(1, L))) then

            write(*, '(I4,3(2x,e15.7))') l, eddi(1, l), xk, xjc(l)

            stop 'elimin eddi(1, l) is nan'

         endif

!        THIS IS AN INGENIOUS RECURSION FORMULA FOR THE SPHERICITY FACTOR
         RLP = RADIUS(L + 1)
         FLP = FL
         FL = three - one / EDDI(1, L)
         RRQ = RRQ * EXP(FL - FLP) * (RL / RLP)**((FLP * RL - FL * RLP) / (RL - RLP))
         EDDI(2, L) = RRQ / RL / RL
         A(1 : JMAX) = WX(1 : JMAX, L)

      ENDDO
     
      CALL MOMENT1(RADIUS(1), NP, P(1 : NP), WX(1 : NP, 1), H)

!     EMFLUX IS THE EMERGENT FLUX, BUT RELATED TO THE INNER RADIUS
      EMFLUX = 4 * H * RADIUS(1) * RADIUS(1)
      EDDI(3, 1) = H / XJC(1)

      deallocate(A)
      deallocate(B)
      deallocate(C)
      deallocate(W)
      deallocate(BX)
      deallocate(WX)

      RETURN

      END SUBROUTINE

      SUBROUTINE SETUP(L, A, B, C, W, JMAX, ND, NP, OPA, ETA, THOMSON, Z, RADIUS, BCORE, DBDR)

!     FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY:
!     SET UP THE MATRICES  A (DIAGONAL), B (FULL), C (DIAGONAL) AND W (VECTOR)

      USE MOMENTS

      implicit none

      integer,intent(in) :: L, ND, NP

      real*8, intent(in) ::                    BCORE, DBDR
      real*8, intent(in), dimension(ND) ::     ETA, OPA, RADIUS, THOMSON
      real*8, intent(in), dimension(ND, NP) :: Z

      integer,intent(out) ::                    JMAX
      real*8, intent(out), dimension(NP) ::     A, C, W
      real*8, intent(out), dimension(NP, NP) :: B

      integer :: J, JMM, JS
      real*8 ::  CORFAC, DT, DTP, DTM, ETAL, G
      real*8 ::  PLUSI, WJG, X, XM, XP, ZLJ, ZLMIN, ZLPLUS

      real*8, parameter :: one = 1.D+0, two = 2.d0, three = 3.D0
     
      JMAX = NP + 1 - L
      JMM =  JMAX - 1
     
C***  EVERY L = 1 ... ND
      X=OPA(L)
      G=-X*THOMSON(L)
      ETAL=ETA(L)
     
C***  MEAN INTENSITY INTEGRATION WEIGHTS FROM SUBROUTINE MOMENT0 (VEKTOR W)
      CALL MOMENT0_SETUP(ND, RADIUS, L, JMAX, Z(1 : ND, 1 : JMAX), W(1 : JMAX))

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

      end subroutine setup

      end module
