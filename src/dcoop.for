      module MOD_DCOOP

      contains

      SUBROUTINE DCOOP (I,DOPA,DETA,XLAMBDA,NF,TL,RNEL,ENTOTL,EN,RSTAR,
     $                  WCHARM,ND,L,NFEDGE,EXPFAC,N,NCHARG,WEIGHT,
     $                  ELEVEL,EION,NOM,EINST,SIGMAKI)

!***********************************************************************
!     DERIVATIVE OF NON-LTE OPACITY AND EMISSIVITY WITH RESPECT TO EN(I)
!     AT CURRENT DEPTH POINT L FOR ALL CONT. FREQUENCIES
!     THE FREQUENCY-DEPENDENT BOUND-FREE CROSS SECTION SIGMA (IN CM**2)
!     IS TAKEN FROM THE ARRAY SIGMAKI
!     IMPORTANT NOTE: THE DERIVATIVES ARE ONLY CALCULATED CORRECTLY FOR
!     THOSE FREQUENCIES AND DEPTH POINTS WHERE THEY ARE NEEDED,
!     I.E. IF (WCHARM(L, K) .GT. 0)
!***********************************************************************

      use MOD_GAUNTFF

      use common_block

      implicit none
      !global variables
      integer, intent(in)                    :: I, L, N, ND
      integer, intent(in)                    :: NF
      real*8,  intent(in)                    :: ENTOTL
      real*8,                    intent(in)  :: RNEL,RSTAR,TL
      integer, dimension(N),     intent(in)  :: NCHARG,NFEDGE,NOM
      real*8,  dimension(N),     intent(in)  :: ELEVEL,EION(N),WEIGHT
      real*8,  dimension(NF),    intent(in)  :: EXPFAC,XLAMBDA
      real*8,  dimension(N,N),   intent(in)  :: EINST
      real*8,  dimension(N),     intent(in)  :: EN
      real*8,  dimension(ND,NF), intent(in)  :: WCHARM
      real*8,  dimension(NF,N),  intent(in)  :: SIGMAKI
      real*8,  dimension(NF),    intent(out) :: DOPA,DETA
	
      !local variables
      real*8, dimension(-1:6)                :: GFF,GIIIX
      real*8                                 :: ABSFAC, C1, C2, C3, CFF, EDGE, EMINDU
      real*8                                 :: EXPEDGE, G, GIII
      real*8                                 :: PRESIG,ROOTTL,SIGMA,SIGMAFF,SUM,T32,TLOG
      real*8                                 :: W,W3,WE,XLAMLOG
      integer                                :: IG, J, K, LOW, NFLOW, NUP
     
!     C1 = H * C / K (CM * ANGSTROEM)
      DATA C1 / 1.4388d0 /
!     C2 = 2 * H * C ( G * CM**3 / S**2)
      DATA C2 / 3.972d-16 /
!     C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      DATA C3 / 2.07d-16 /
!     CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION (ALLEN PAGE 100)
      DATA CFF / 1.370d-23 /
     
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL
     
      DETA(:) = .0d0
      DOPA(:) = .0d0

C*** IF DERIVATIVE WITH RESPECT TO ELECTRON DENSITY REQUIRED GOTO SPECIAL BRANCH
      IF (I .EQ. N+1) GOTO 3
     
C***  BOUND-FREE TRANSITIONS INVOLVING LEVEL I  ********************************
      DO 2, J=1,N
      	IF (I .EQ. J) goto 2
C***  LEVELS MUST BELONG TO THE SAME ELEMENT
      IF (NOM(I) .NE. NOM(J)) GOTO 2
     
      IF (I.LT.J) THEN
            LOW=I
            NUP=J
            ELSE
            LOW=J
            NUP=I
      ENDIF
C***  UPPER LEVEL MUST BE GROUND STATE OF THAT ION
      IF (NCHARG(NUP) .EQ. NCHARG(NUP-1)) GOTO 2
C***  CHARGES MUST DIFFER BY 1
      IF (NCHARG(NUP).NE.NCHARG(LOW)+1 ) GOTO 2
C***  SEARCH FOR THE INDEX OF IONIZATION EDGE
      EDGE=EION(LOW)-ELEVEL(LOW)
      EXPEDGE=EXP(C1*EDGE/TL)
      NFLOW=NFEDGE(LOW)
     
      IF (I .LT. J) THEN
C***  I IS LOWER LEVEL
      	DO K = 1, NFLOW

      		DOPA(K) = DOPA(K) + SIGMAKI(K, idx_nlte(LOW))

      	ENDDO
     
      ELSE
C***  I IS UPPER LEVEL
      	WE=C3*RNEL*ENTOTL/T32 *WEIGHT(LOW)/WEIGHT(NUP)*EXPEDGE
      	DO K=1,NFLOW
	      	W=1.d8/XLAMBDA(K)
C*** RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
            G=WE*EXPFAC(K)
            SIGMA=SIGMAKI(K, idx_nlte(LOW))
            DOPA(K)=DOPA(K)-G*SIGMA
            DETA(K)=DETA(K)+G*SIGMA
        ENDDO
      ENDIF
     
 2    ENDDO
     
     
C***  FREE-FREE CONTRIBUTION OF LEVEL I  ***************************************
      IF (NCHARG(I) .LE. 0) GOTO 4
      TLOG=LOG10(TL)
      DO 20 K=1,NF
      IF (WCHARM(L,K) .EQ. .0 ) GOTO 20
		W=1.d8/XLAMBDA(K)
		W3=W*W*W
		XLAMLOG=LOG10(XLAMBDA(K))
		SIGMAFF=FLOAT(NCHARG(I)*NCHARG(I))*CFF/W3/ROOTTL
C***  FREE-FREE GAUNT FACTOR GIII CALCULATED
		CALL GFFLOG (GIII,NCHARG(I),XLAMLOG,TLOG)
		SIGMAFF=SIGMAFF*GIII
		SUM=RNEL*ENTOTL*SIGMAFF
		EMINDU=SUM*EXPFAC(K)
		SUM=SUM-EMINDU
		DOPA(K)=DOPA(K)+SUM
		DETA(K)=DETA(K)+EMINDU
   20 CONTINUE
      GOTO 4
     
C***  I .EQ. N+1 : DERIVATIVE WITH RESPECT TO ELECTRON DENSITY  ================
    3 CONTINUE
C***  BOUND-FREE OPACITIES (DERIVATIVE TO ELECTRON DENSITY) ********************
      DO 7 NUP=2,N
C***  UPPER LEVEL MUST BE GROUND STATE OF HIGHER IONIZATION STAGES
      	IF ((NOM(NUP) .NE. NOM(NUP-1)) .OR.
     $    (NCHARG(NUP) .NE. NCHARG(NUP-1))) GOTO 7
      DO 5 LOW=1,NUP-1
C***  LEVELS MUST BELONG TO THE SAME ELEMENT
      IF (NOM(LOW) .NE. NOM(NUP)) GOTO 5
C***  CHARGES MUST DIFFER BY 1
      IF (NCHARG(NUP).NE.NCHARG(LOW)+1 ) GOTO 5
          EDGE=EION(LOW)-ELEVEL(LOW)
          NFLOW=NFEDGE(LOW)
          EXPEDGE=EXP(C1*EDGE/TL)
          !***  NOTE: FACTOR RNEL OMITTED (DERIVATIVE TO ELECTRON DENSITY)
          WE=C3*     ENTOTL/T32 *WEIGHT(LOW)/WEIGHT(NUP)*EXPEDGE
          DO K=1,NFLOW
            W=1.d8/XLAMBDA(K)
            SIGMA=SIGMAKI(K, idx_nlte(LOW))
            !G DIVIDED BY RNEL
            G=WE*EXPFAC(K)
            DETA(K)=DETA(K)+G*SIGMA*EN(NUP)
            DOPA(K)=DOPA(K)-G*SIGMA*EN(NUP)
         ENDDO
     
    5 CONTINUE
    7 CONTINUE
     
C***  FREE-FREE OPACITIES (DERIVATIVE TO ELECTRON DENSITY) *********************
C***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      TLOG=LOG10(TL)
      DO K=1,NF
        IF (WCHARM(L,K) .ne. .0d0 ) then
          W=1.d8/XLAMBDA(K)
          W3=W*W*W
          XLAMLOG=LOG10(XLAMBDA(K))
          DO IG = 0,6
            IF (IG.GT.0) THEN
              CALL GFFLOG (GIIIX(IG),IG,XLAMLOG,TLOG)
            ELSE
              GIIIX(IG)=.0d0
            ENDIF
            GFF(IG)=GIIIX(IG)*IG*IG
          ENDDO
C***      PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
          PRESIG =CFF/W3/ROOTTL
          DO J=1,N
            SIGMAFF   = PRESIG*GFF(NCHARG(J))
            SUM       = ENTOTL*EN(J)*SIGMAFF
            EMINDU    = SUM*EXPFAC(K)
            SUM       = SUM-EMINDU
            DOPA(K)   = DOPA(K)+SUM
            DETA(K)   = DETA(K)+EMINDU
          ENDDO
        endif
      ENDDO
     
    4 CONTINUE
      DO K=1,NF
        W=1.d8/XLAMBDA(K)
C       IF (XLAM.GT.227.83 .AND. XLAM.LE.3100 ) THEN
C         ELDEN=RNEL*ENTOTL
C         CALL LINSCA (XLAM,ELDEN,SCAFAC,ABSFAC)
C***        MULTIPLY THE TRUE ABSORPTION TO ACCOUNT FOR LINE ABSORPTION
C***        MULTIPLY SUM TO ACOUNT FOR LINE-SCATTERING
C       ELSE
        ABSFAC=1.d0
C       ENDIF
        W3=W*W*W
        DOPA(K)=DOPA(K)*ENTOTL*RSTAR*ABSFAC
        DETA(K)=DETA(K)*ENTOTL*RSTAR*C2*W3
      ENDDO
     
      RETURN
      END subroutine
      end module
