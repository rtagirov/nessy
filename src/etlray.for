      MODULE MOD_ETLRAY

      CONTAINS

      SUBROUTINE ETLRAY(U,Z,OPA,OPAL,ETA,ETAL,XJLMEAN,RADIUS,ND,NP,JP,P,
     $                  UB,GA,H,QQ,S,V,VA,VB,VELO,GRADI,PP,
     $                  BMHO,BMNO,BMHI,BMNI,BCORE,DBDR,
     $                  XJCIND,ELEVEL,LOW,
     $                  PHI,PWEIGHT,NFL,DELTAX,XJ,XH,XK,XN,W0,W1,W2,W3,
     $                  NL, XLAM, LO)

!     LINE RADIATION TRANSFER IN THE COMOVING FRAME FROM A GIVEN SOURCE FUNCTION
!     FOR A GIVEN IMPACT-PARAMETER, THE INTEGRATION IS CARRIED OUT IN SPACE
!     AND FREQUENCY AND RESULTS ARE INTEGRATED UP TO YIELD THE MEAN INTENSITY
!     OF THE LINE, XJLMEAN = J-NUE-BAR, AND THE 0 TO 3 MOMENTS J, H, K, N

!     RINAT TAGIROV:
!     See Fig. 7-29 in Mihalas, Stellar Atmospheres, 2nd edition, 1978
!     for the geometry of the ray-by-ray solution as well as
!     Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows. 
!     I.Computational method for equivalent-two-level-atom source functions" (MKH)

      use MOD_MDV
      use MOD_MOMADD
      use MOD_COFREQ
      use MOD_CMFSET
      use MOD_VMALV
      use MOD_INVTRI
      use MOD_GMALU
      USE MOD_CALCLAMBDAS

      USE FILE_OPERATIONS

      use matoper

      IMPLICIT REAL*8(A - H, O - Z)

      parameter (one = 1.0D0)
     
      DIMENSION U(ND),XJLMEAN(ND),RADIUS(ND),V(ND),Z(ND,NP)
      DIMENSION VELO(ND),GRADI(ND),PP(ND)
      DIMENSION PHI(NFL),PWEIGHT(NFL)
      DIMENSION BMHO(NFL),BMNO(NFL),BMHI(NFL),BMNI(NFL)
      DIMENSION XJ(NFL,ND),XH(NFL,ND),XK(NFL,ND),XN(NFL,ND)
      DIMENSION W0(ND),W1(ND),W2(ND),W3(ND)
      DIMENSION OPA(ND),OPAL(ND),ETA(ND),ETAL(ND)
      DIMENSION XJCIND(ND),ELEVEL(2)
      DIMENSION H(*), UB(*),S(*),QQ(*), P(*)

      DIMENSION VA(ND), VB(ND), GA(ND)

!     RINAT TAGIROV:
!     TA, TB, TC: see CMFSET.FOR
!     TAL, TBL, TCL are TA, TB and TC, but without the contribution from continuum.
!     They were introduced to implement the local approximate Lambda-operator for lines.
!     OPA = ZERO, ETA = ZERO for TAL, TBL, TCL (first call of CMFSET subroutine).
!     If OPA and ETA are not equal to zero it renders the calculation of the local approximate Lambda-operator impossible.

      REAL*8, ALLOCATABLE, DIMENSION(:) :: TA, TB, TC
      REAL*8, ALLOCATABLE, DIMENSION(:) :: TAL,TBL,TCL
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ZERO

      logical THIN

!     Diagonal of the inverse of TRI(TAL, TBL, TCL), which is matrix T in MKH formulation (equation 2.39), for a particular frequency within a line
      REAL*8, ALLOCATABLE, DIMENSION(:) :: INV_T_DIAG

!     INV_T_DIAG averaged over the line profile and over all angles, i.e. the Local approximate lambda-Operator for a given line
      REAL*8, DIMENSION(ND), INTENT(INOUT) :: LO

      COMMON / COMFUN / DELTAV,XMIN

      LMAX = MIN0(NP + 1 - JP, ND)
      LZ = LMAX - 1

      ALLOCATE(ZERO(ND))

      ALLOCATE(INV_T_DIAG(LMAX))

      ALLOCATE(TA(LMAX))
      ALLOCATE(TB(LMAX))
      ALLOCATE(TC(LMAX))

      ALLOCATE(TAL(LMAX))
      ALLOCATE(TBL(LMAX))
      ALLOCATE(TCL(LMAX))

      ZERO(1 : ND) = 0.0D0
      INV_T_DIAG(1 : LMAX) = 0.0D0

      TA(1 : LMAX) = 0.0D0
      TB(1 : LMAX) = 0.0D0
      TC(1 : LMAX) = 0.0D0

      TAL(1 : LMAX) = 0.0D0
      TBL(1 : LMAX) = 0.0D0
      TCL(1 : LMAX) = 0.0D0

C***  VERSION OF OUTER BOUNDARY CONDITION
C***  0 : NO INCIDENT RADIATION
C***  1 : IMINUS = S  WITHIN THE SCHARMER CORES
C***  2 : IMINUS = S * (one-EXP(-TBOUND))
C***  3 : IMINUS = S * (one-EXP(-TBOUND)), S IN SOBOLEV APPROXIMATION
      IVERSION = 0
C***  NON-ZERO INCIDENT RADIATION RESTRICTED TO RESONANCE LINES (OPTIONALLY)
CCC   IF (ELEVEL(LOW) .NE. .0) IVERSION=0
     
C***  BLUE-WING BOUNDARY CONDITION..  U AND V SPECIFIED
C***  THE FEAUTRIER-FLUX V (AT DEPTH INTERSTICES) IS DERIVED FROM U VIA THE DGL.
      DO 5 L=1,LZ
      X=0.5d0*(OPA(L)+OPA(L+1))
    5 V(L)=(U(L+1)-U(L))/X/(Z(L,JP)-Z(L+1,JP))
     
C***  PP(L) = VELOCITY GRADIENT, PROJECTED ON THE PRESENT RAY J, /DELTAX
      DO 4 L=1,LMAX
      RL=RADIUS(L)
      Y=Z(L,JP)/RL
      YY=Y*Y
    4 PP(L)=(YY*GRADI(L)+(one-YY)*VELO(L)/RL)/DELTAX
     
C***  GENERATING INTEGRATION WEIGHTS APPROPRIATE FOR THE CONSIDERED RAY
C***  WHICH WILL BE USED TO INTEGRATE THE 0. TO 3. MOMENTS
      CALL       MOMADD (JP,LMAX,ND,NP,Z,RADIUS,P,W0,W1,W2,W3,
     $                   WBMHO,WBMNO,WBMHI,WBMNI)
     
      IF (IVERSION .EQ. 1) THEN
C***  MODIFIED OUTER BOUNDARY CONDITION:
C*** IN THE LINE CORES, IMINUS=S(1) IS ASSUMED
C***  THE LINE CORE IS DEFINED AS IN THNE SCHARMER PAPERS WITH GAMMA=1.
         DELTAV=VELO(1)-VELO(ND)
         XRED=.0d0
         XBLUE=.0d0
         XMAX=DELTAX*FLOAT(NFL-1)*.5d0
         IF (OPAL(1) .LE. .0d0) GOTO 3
         GDT=PP(1)*DELTAX/OPAL(1)
         ERXMIN=ERF_INF(-XMAX)
         CALL COFREQ (XR, XB, XMAX, ERXMIN, GDT, THIN)
         XRED=-XB
         XBLUE=-XR
    3    CONTINUE
      ENDIF
     
C***  LOOP FOR ALL (ORIGINAL) FREQUENCY POINTS
C***  THE CURRENT INDEX K REFERS TO THE RESULTING U,V , WHERE ALL QUANTITIES
C***  ARE TAKEN. THE FREQUENCY DERIVATIVE STRETCHES BACK TO K-1
     
      DO 1 K = 1, NFL
 
C***  OUTER BOUNDARY CONDITIONS
      IF (IVERSION .EQ. 0) XIMINUS=0.
      IF (IVERSION .EQ. 1) THEN
      XL=XMAX-(K-1)*DELTAX
      IF (XL .GT. XRED .AND. XL .LT. XBLUE) THEN
            AK=OPA(1)+PHI(K)*OPAL(1)
            EK=ETA(1)+PHI(K)*ETAL(1)
            XIMINUS=EK/AK
            ELSE
            XIMINUS=.0d0
            ENDIF
            ENDIF
      IF (IVERSION .EQ. 2) THEN
            AK=OPA(1)+PHI(K)*OPAL(1)
            EK=ETA(1)+PHI(K)*ETAL(1)
            TBOUND=RADIUS(1)*AK
            XIMINUS=(one-EXP(-TBOUND))*EK/AK
            ENDIF
      IF (IVERSION .EQ. 3) THEN
            GR=GRADI(1)
            GT=VELO(1)/RADIUS(1)
            SSOBO=3.*XJCIND(1)*(GR/GT)*
     *                  (one-EXP(-OPAL(1)/GR))/(one-EXP(-OPAL(1)/GT))
            ETASOBO=OPAL(1)*SSOBO
            EK=ETA(1)+PHI(K)*ETASOBO
            AK=OPA(1)+PHI(K)*OPAL(1)
            TBOUND=RADIUS(1)*AK
            XIMINUS=(one-EXP(-TBOUND))*EK/AK
            ENDIF
     
C***  FREQUENCY DERIVATIVE OF IMINUS
      IF (K .GT. 1) DXI=XILASTK-XIMINUS
      XILASTK=XIMINUS
     
      IF (K .EQ. 1) GOTO 2

      CALL CMFSET(PHI(K),Z(1,JP),
     $            ND,LMAX,TAL,TBL,TCL,UB,VA,VB,GA,H,S,ZERO,OPAL,ZERO,
     $            ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

      CALL CMFSET(PHI(K),Z(1,JP),
     $            ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,OPAL,ETA,
     $            ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

      INV_T_DIAG = INVTRIDIAG(TAL, TBL, TCL)

!      do iii = 1, LMAX

!         write(*, '(A,4(5x,I4),5x,e15.7)'), 'etlray:', NL, JP, K, iii, inv_t_diag(iii)
!         write(*, '(A,3(5x,I4),5x,e15.7,5x,i4,5x,e15.7)'), 'etlray:', NL, JP, K, pweight(k), iii, w0(iii)

!      enddo

!      IF (NL .EQ. 1) THEN
!      IF ((K .EQ. 31) .AND. ((JP .EQ. 1) .OR. (JP .EQ. 54))) THEN

!         CALL OPEN_TO_APPEND(14, 'line_feautrier_matrix.out')

!         DO L = 1, LMAX

!            WRITE(14, 10) NL, JP, K, L, TAL(L), TBL(L), TCL(L), INV_T_DIAG(L)

!         ENDDO

!      ENDIF

      CALL VMALV (VA,VB,V,QQ,LMAX)
      CALL VADD (QQ,S,LMAX) ! QQ = QQ + S
      CALL MDV (UB,U,LMAX)
      CALL VADD (U,QQ,LMAX) ! U = U + QQ
      CALL INVTRI (TA,TB,TC,U,LMAX)
C***  NOW U IS THE FIELD AT THE NEW INDEX K
      CALL MDV (H,V,LZ)
      CALL GMALU (GA,U,S,LMAX)
      CALL VADD (V,S,LZ) ! V = V + S

C***  ADDING THE NEW U(L) TO THE MEAN INTENSITY XJMEAN
C***  AND TO THE 0. AND 2. MOMENTS XJ AND XK
    2 DO 8 L = 1, LMAX

      IF (INV_T_DIAG(L) .NE. INV_T_DIAG(L))
     $ WRITE(*, '(A,1x,I3,1x,I3,1x,I3,1x,E23.15)')
     $ 'ETLRAY: INV_TRI_DIAG IS NaN',
     $ K, JP, L, INV_T_DIAG(L)

      IF (INV_T_DIAG(L) .LT. 0.0D0)
     $ WRITE(*, '(A,1x,I3,1x,I3,1x,I3,1x,E23.15)')
     $ 'ETLRAY: INV_TRI_DIAG IS NEGATIVE',
     $ K, JP, L, INV_T_DIAG(L)

      LO(L) = LO(L) + INV_T_DIAG(L) * W0(L) * PWEIGHT(K)

      XJLMEAN(L) = XJLMEAN(L) + U(L) * W0(L) * PWEIGHT(K)

      XJ(K,L)=XJ(K,L)+W0(L)*U(L)
    8 XK(K,L)=XK(K,L)+W2(L)*U(L)
      
C***  ADDING THE NEW V TO THE 1. AND 3. MOMENTS XH AND XN
      DO 9 L=1,LZ
      XH(K,L)=XH(K,L)+W1(L)*V(L)
    9 XN(K,L)=XN(K,L)+W3(L)*V(L)
C***  ADDING U(1) TO THE BOUNDARY MOMENTS H AND N (BMH, BMN)
C***  MODIFICATION FOR NONZERO INCIDENT RADIATION FROM TRUNCATED LAYERS
      BMHO(K)=BMHO(K)+WBMHO*(U(1)-XIMINUS)
      BMNO(K)=BMNO(K)+WBMNO*(U(1)-XIMINUS)
      IF (LMAX.LT.ND) GOTO 1
      BMHI(K)=BMHI(K)+WBMHI*U(ND)
      BMNI(K)=BMNI(K)+WBMNI*U(ND)
    1 CONTINUE

      DEALLOCATE(ZERO)

      DEALLOCATE(INV_T_DIAG)

      DEALLOCATE(TA)
      DEALLOCATE(TB)
      DEALLOCATE(TC)

      DEALLOCATE(TAL)
      DEALLOCATE(TBL)
      DEALLOCATE(TCL)

      RETURN

10    FORMAT(I3,8x,I2,8x,I2,8x,I2,8x,E15.7,8x,E15.7,8x,E15.7,8x,E15.7)

      END SUBROUTINE

      END MODULE
