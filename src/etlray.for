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

      use MOD_MOMADD
      use MOD_COFREQ
      use MOD_INVTRI
      use MOD_GMALU
      USE MOD_CALCLAMBDAS
      USE FILE_OPERATIONS
      USE UTILS
      USE MATOPER

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
      CALL MOMADD(JP,LMAX,ND,NP,Z,RADIUS,P,W0,W1,W2,W3,WBMHO,WBMNO,WBMHI,WBMNI)
     
C***  LOOP FOR ALL (ORIGINAL) FREQUENCY POINTS
C***  THE CURRENT INDEX K REFERS TO THE RESULTING U,V , WHERE ALL QUANTITIES
C***  ARE TAKEN. THE FREQUENCY DERIVATIVE STRETCHES BACK TO K-1
     
      DO 1 K = 1, NFL
 
C***  OUTER BOUNDARY CONDITIONS
      XIMINUS=0.0d0
     
C***  FREQUENCY DERIVATIVE OF IMINUS
      IF (K .GT. 1) DXI = XILASTK - XIMINUS

      XILASTK = XIMINUS
     
      IF (K .EQ. 1) GOTO 2

      CALL CMFSET(PHI(K),Z(1,JP),ND,LMAX,TAL,TBL,TCL,UB,VA,VB,GA,H,S,ZERO,OPAL,ZERO,
     $            ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

      INV_T_DIAG = INVTRIDIAG(TAL, TBL, TCL)

      call assert(.not. any(isnan(inv_t_diag)),     'etlray: inv_t_diag is NaN')
      call assert(.not. any(inv_t_diag .lt. 0.0d0), 'etlray: inv_t_diag is negative')

      CALL CMFSET(PHI(K),Z(1,JP),ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,OPAL,ETA,
     $            ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

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

      SUBROUTINE CMFSET(PHIK,Z,ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,OPAL,
     $                  ETA,ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

!     THIS SUBROUTINE IS TO SET UP THE ARRAY ELEMENTS FOR THE CMF FORMALISM

!     RINAT TAGIROV:
!     It follows formalism described in Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows. 
!     I.Computational method for equivalent-two-level-atom source functions" (MKH),
!     but with the source function of general form, i.e the following equation is solved:
!     T_{k,j}u_{k,j} + U_{k,j}u_{k-1,j} + V_{k,j}v_{k-1,j} = S_{k,j}
!     instead of Eq. (2.39) in MKH.
!     Specifically see their Eq. (2.31) and Appendix C.
!     TA, TB, TC are the diagonals of matrix T_{k,j} (tridiagonal).
!     VA, VB are the diagonals of matrix V_{k,j} (bidiagonal).
!     UB is the matrix U_{k,j} (diagonal).
!     H is the matrix H_{k,j} (diagonal).
!     GA is matrix G_{k,j} which is bidiagonal but its diagonal element 
!     is minus one multiplied by the non-diagonal one, so only GA array is needed here.
!     The signs do not match their formulation because apparently the routine is written
!     for the aforementioned matrix equation having slightly different 
!     form with some signs taken out of the matrix elements
!     and put directly into the equation.

      implicit real*8(a-h,o-z)

      parameter (half=0.5d0, one=1.0d0, two=2.0d0)

      DIMENSION S(ND),OPA(ND),OPAL(ND),ETA(ND),ETAL(ND),VA(ND),VB(ND)
      DIMENSION UB(ND),GA(ND),H(ND)
      DIMENSION PP(ND),Z(ND)

      DIMENSION TA(LMAX), TB(LMAX), TC(LMAX)

      LZ = LMAX - 1
     
!     OUTER BOUNDARY CONDITION - FIRST ORDER
      AK=OPA(1)+PHIK*OPAL(1)
      AZ=half*(AK+OPA(2)+PHIK*OPAL(2))
      TAUZ=Z(1)-Z(2)
      DX=PP(1)/AK
      S(1)=XIMINUS+DX*DXI
      TC(1)=one/AK/TAUZ
      TB(1)=TC(1)+DX+one
      UB(1)=DX
      VB(1)=0.0

!     FOR G AND H THE MATRIX ELEMENT S ARE NOT DIFFERENT FROM INNER POINTS
      DTZM=one/AZ/TAUZ
      DXZM=(PP(1)+PP(2))/AZ/two
      DAZM=DTZM/(one+DXZM)
      DBZM=DXZM/(one+DXZM)
      GA(1)=-DAZM
      H(1)=DBZM
      IF(LZ.LT.2) GOTO 2
     
!     NON - BOUNDARY POINTS
      DO 1 L=2,LZ
      AK=OPA(L)+PHIK*OPAL(L)
      EK=ETA(L)+PHIK*ETAL(L)
      S(L)=EK/AK
      AZ=half*(AK+OPA(L+1)+PHIK*OPAL(L+1))
      AZM=half*(AK+OPA(L-1)+PHIK*OPAL(L-1))
      TAU=half*(Z(L-1)-Z(L+1))
      TAUZ=Z(L)-Z(L+1)
      TAUZM=Z(L-1)-Z(L)
      DT=one/AK/TAU
      DTZ=one/(AZ*TAUZ)
      DTZM=one/(AZM*TAUZM)
      DX=PP(L)/AK
      DXZ=(PP(L)+PP(L+1))*half/AZ
      DXZM=(PP(L)+PP(L-1))*half/AZM
      DAZ=DTZ/(one+DXZ)
      DAZM=DTZM/(one+DXZM)
      DBZ=DXZ/(one+DXZ)
      DBZM=DXZM/(one+DXZM)

      TA(L)=DT*DAZM
      TC(L)=DT*DAZ
      TB(L)=TA(L)+TC(L)+DX+one

      UB(L)=DX

      VA(L)=-DT*DBZM
      VB(L)=DT*DBZ

      GA(L)=-DAZ

      H(L)=DBZ

    1 CONTINUE
     
    2 L = LMAX
      IF (LMAX .LT. ND) GOTO 4
     
!     INNER BOUNDARY CONDITION (CORE RAYS) - ONLY TO FIRST ORDER
      AK=OPA(ND)+PHIK*OPAL(ND)
      S(ND)=BCORE+DBDR*Z(ND)/AK
      TAUZ=Z(ND-1)-Z(ND)
      DT=one/TAUZ/AK
      DX=PP(L)/AK
      TA(L)=DT
      TB(L)=DT+DX+one
      UB(L)=DX
      VA(L)=0.0d0

      RETURN
     
!     INNER BOUNDARY CONDITION (NON-CORE RAYS) - SECOND ORDER
    4 AK=OPA(LMAX)+PHIK*OPAL(LMAX)
      EK=ETA(L)+PHIK*ETAL(L)
      S(L)=EK/AK
      TAUZ=Z(LZ)
      DT=one/AK/TAUZ
      DX=PP(L)/AK
      DA=DT/(one+DX)
      DB=DX/(one+DX)
      TA(L)=two*DT*DA
      TB(L)=TA(L)+DX+one
      UB(L)=DX
      VA(L)=-two*DT*DB

      RETURN

      END SUBROUTINE

      END MODULE
