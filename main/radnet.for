      MODULE MOD_RADNET

      CONTAINS

!**********  MODULNAME: RADNET    ******* 06/08/87  20.15.04.******   109 KARTEN
      SUBROUTINE RADNET(NDIM,N,ENLTE,TL,WEIGHT,NCHARG,EION,ELEVEL,EINST,
     $                  SL,EN,NOM,RRATE,XLAMBDA,FWEIGHT,
     $                  XJC,NF,L,XJL,ND,SIGMAKI,LASTIND,
     $                  LEVEL, DP, JOBNUM, ITNEL)

      !*****************************************************************
      !***  RADIATIVE RATE COEFFICIENT MATRIX RRATE IS CALCULATED AT
      !***  DEPTH POINT L FROM THE GIVEN RADIATION FIELD.
      !***  NOTE THAT THIS SUBROUTINE CALCULATES NETTO RATE COEFFICIENTS
      !***  FOR NON-RUDIMENTAL LINE TRANSITIONS
      !*****************************************************************

      !***  NOTE: ND = 1 and L = 1 ! See the call of radnet.for from coma.for

      USE MOD_ISRCHFGT
      USE MOD_XRUDI
      USE UTILS
      USE PHYS
      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT REAL * 8(A - H, O - Z)

      PARAMETER (one = 1.0D0)

      DIMENSION EINST(NDIM,NDIM)
      DIMENSION XJC(ND,NF),XJL(ND,LASTIND)
      DIMENSION NCHARG(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),ENLTE(NDIM),WEIGHT(NDIM),EN(NDIM)
      DIMENSION NOM(N)
      DIMENSION SIGMAKI (NF,N)
      DIMENSION XLAMBDA(NF), FWEIGHT(NF)

      INTEGER, INTENT(IN) :: DP

      REAL*8 :: XJ_LTE, EMINDU_LTE

      REAL*8, DIMENSION(NDIM, NDIM), INTENT(OUT) :: RRATE

      INTEGER, INTENT(IN) :: JOBNUM, ITNEL

      REAL*8, DIMENSION(LASTIND), INTENT(IN) :: SL

      CHARACTER*10, DIMENSION(NDIM), INTENT(IN) :: LEVEL

      REAL*8 :: FREQ_FACT

      REAL*8 :: XJCLK_LTE

      LOGICAL :: ARR_LTE_COND, DEPTH_COND

      REAL*4 ::  JSRATIO

      !***  C1 = H * C / K      (CM * KELVIN)
      DATA C1 / 1.4388d0 /
      !***  C2 = 2 * H * C      (H AND C IN CGS UNITS)
      DATA C2 / 3.9724d-16 /
      !***  C3 = 4 * PI / H / C (CGS UNITS)
      DATA C3 / 0.06327d18 /
      !***  PI8 = 8 * PI
      DATA PI8 / 25.13274123d0 /

      DEPTH_COND = (DP .LE. 91 .AND. DP .GE. 80) .OR. (DP .LE. 11)

      ARR_LTE_COND = ALLOCATED(ARR_LTE) .AND. ITNEL .EQ. 1

      call assert(ND*L==1, 'ND and L must both be 1 '//int2str(ND)//int2str(L))

      ARR(DP, 1 : N, 1 : N) = 0.0D0
      RBR(DP, 1 : N, 1 : N) = 0.0D0

      IF (ARR_LTE_COND) ARR_LTE(DP, 1 : N, 1 : N) = 0.0D0

      !***  LOOP OVER ALL TRANSITIONS  ---------------------------------

      IND = 0

      DO NUP = 2, N

        LLOW : DO LOW = 1, NUP - 1

          IF (NOM(LOW) .NE. NOM(NUP)) GOTO 14

          IF (NCHARG(LOW) .NE. NCHARG(NUP)) GOTO 8

          !***  LINE TRANSITION   **************************************

          IND = IND + 1

          WAVENUM = ELEVEL(NUP) - ELEVEL(LOW)
          W3 = WAVENUM * WAVENUM * WAVENUM

          !***  CHECK WHETHER THIS TRANSITION IS ONLY RUDIMENTAL

!***************************************************************************************************************************************
!***************************************************************************************************************************************
! THAT IS HOW IT'S DONE WITH THE RADIATIVE BRACKET

          IF (EINST(LOW, NUP) .NE. -2.0D0 .AND. SL(IND) .NE. 0.0D0) THEN

             XJ = XJL(L, IND)

             JSRATIO = XJ / SL(IND)

!            RINAT TAGIROV:
!            At the innermost point the value of radiative bracket seems to be wrong. It abruptly increases.
!            For some reason the XJ / SL(IND) ratio there is much further from one than at the preceding point.
!            For possible reasons see the commentary in ETL.FOR about the boundary conditions and weights.
!            The explanation given there however does not seem to be in line with the fact that the value
!            of radiative bracket at the outermost point looks correct.
!            Though this may be attributed to insensitivity of the radiative bracket to the outer boundary condition.
!            But then again I can not think of any reason for such insensitivity.
             RRATE(LOW, NUP) = 0.0D0
             RRATE(NUP, LOW) = EINST(NUP, LOW) * (one - XJ / SL(IND))

!             RRATE(NUP, LOW) = EINST(NUP, LOW) * (one - JSRATIO)

!             IF (IND .EQ. 1) RRATE(NUP, LOW) = EINST(NUP, LOW) * (one - XJ / SL(IND)) + 1000 * 0.25 * 8.2249

          ELSEIF (EINST(LOW, NUP) .NE. -2.0D0 .AND. SL(IND) .EQ. 0.0D0) THEN

             XJ = XJL(L, IND)

             EMINDU = EINST(NUP, LOW) * XJ / C2 / W3

             RRATE(LOW, NUP) = EMINDU * WEIGHT(NUP) / WEIGHT(LOW)

             RRATE(NUP, LOW) = EINST(NUP, LOW) + EMINDU

          ELSE

          !*** TRANSITION IS RUDIMENTAL -- RADIATION FIELD FROM INTERPOLATION OF CONT.

             CALL XRUDI(XJ, WAVENUM, XJC, XLAMBDA, ND, NF, L)

             EMINDU = EINST(NUP, LOW) * XJ / C2 / W3

             RRATE(LOW, NUP) = EMINDU * WEIGHT(NUP) / WEIGHT(LOW)

             RRATE(NUP, LOW) = EINST(NUP, LOW) + EMINDU

          ENDIF

!**************************************************************************************************************************************
!**************************************************************************************************************************************

!***************************************************************************************************************************************
!***************************************************************************************************************************************
! THAT IS HOW IT'S DONE WITHOUT THE RADIATIVE BRACKET (TO ACTUALLY REMOVE THE RADIATIVE BRACKET FROM THE CALCULATION REPLACE 
! ARR ARRAY HERE WITH THE RRATE ARRAY AND COMMENT OUT THE IMPLEMENTATION OF THE RADIATIVE BRACKET ABOVE)

          IF (EINST(LOW, NUP) .NE. -2.0D0) THEN

             XJ = XJL(L, IND)

	     IF (ARR_LTE_COND) XJ_LTE = XJL_LTE(DP, IND)

          ELSE

          !*** TRANSITION IS RUDIMENTAL -- RADIATION FIELD FROM INTERPOLATION OF CONT.

             CALL XRUDI(XJ, WAVENUM, XJC, XLAMBDA, ND, NF, L)

             IF (ARR_LTE_COND) CALL XRUDI(XJ_LTE, WAVENUM, XJC_LTE(DP, 1 : NF), XLAMBDA, ND, NF, L)

          ENDIF

          EMINDU = EINST(NUP, LOW) * XJ / C2 / W3

          ARR(DP, LOW, NUP) = EMINDU * WEIGHT(NUP) / WEIGHT(LOW)

          ARR(DP, NUP, LOW) = EINST(NUP, LOW) + EMINDU

          IF (ARR_LTE_COND) XJ_LTE = PLANCK_FUNC(WAVENUM * light_speed, TL) ! RINAT TAGIROV: ONE HAS TO CALCULATE THE LTE ABSOLUTE RADIATIVE RATES 
                                                                            ! WITH THE PLANCK FUNCTION TO COMPARE THEM WITH THE NLTE ONES

	  IF (ARR_LTE_COND) EMINDU_LTE = EINST(NUP, LOW) * XJ_LTE / C2 / W3
	  IF (ARR_LTE_COND) ARR_LTE(DP, LOW, NUP) = EMINDU_LTE * WEIGHT(NUP) / WEIGHT(LOW)
	  IF (ARR_LTE_COND) ARR_LTE(DP, NUP, LOW) = EINST(NUP, LOW) + EMINDU_LTE

!**************************************************************************************************************************************
!**************************************************************************************************************************************

          CYCLE LLOW

  8       CONTINUE
          !***  CHARGE DIFFERENCE MUST BE 1
          IF (NCHARG(NUP) .NE. NCHARG(LOW)+1 ) GOTO 14
          !***  UPPER LEVEL MUST BE A GROUND STATE
          IF (NCHARG(NUP) .NE. NCHARG(NUP-1)+1) GOTO 14

          !***  CONTINUUM TRANSITION (NET RADIATIVE BRACKETS) ****************
          !***  SIGMAKI = PRECALCULATED CROSS SECTION IN CM**2
          !***  EDGE = THRESHOLD ENERGY IN KAYSER *******
          EDGE = EION(LOW) - ELEVEL(LOW)
          EDGELAM = 1.0D8 / EDGE

          !***  RATE INTEGRAL
          REC = 0.0D0
          !***  FIND EDGE FREQUENCY INDEX
          NFEDGE = ISRCHFGT(NF, XLAMBDA, 1, EDGELAM) - 1

          L2 : DO K = 1, NFEDGE

!           CALCULATION OF THE RADIATIVE BRACKET

            WAVENUM = 1.0D8 / XLAMBDA(K)
            W2 = WAVENUM * WAVENUM
            W3 = W2 * WAVENUM

            XJCLK = XJC(L, K)

            SIGMA = SIGMAKI(K, LOW)

            !***  CALCULATE BOUND-FREE SOURCE FUNCTION FOR TRANSITION LOW-UP ONLY
            EXFAC = EXP(-C1 * WAVENUM / TL)
            G = EXFAC * ENLTE(LOW) / ENLTE(NUP)
            SBF = C2 * W3 / (EN(LOW) / (EN(NUP) * G) - one)
            REC = REC + SIGMA * W2 * EXFAC * (one - XJCLK / SBF) * FWEIGHT(K)

!           CALCULATION OF THE ABSOLUTE RADIATIVE RATES (AS GIVEN BY RUTTEN "Radiative Transfer In Stellar Atmospheres")

            FREQ_FACT = (SIGMA / WAVENUM) * FWEIGHT(K)

!           BOUND - FREE:
            ARR(DP, LOW, NUP) = ARR(DP, LOW, NUP) + FREQ_FACT * XJCLK ! FWEIGHT is the frequency interval (according to the trapezoidal integration rule)

!            IF (ARR_LTE_COND) XJCLK_LTE = XJC_LTE(DP, K)
            IF (ARR_LTE_COND) XJCLK_LTE = PLANCK_FUNC(WAVENUM * light_speed, TL)

            IF (ARR_LTE_COND) ARR_LTE(DP, LOW, NUP) = ARR_LTE(DP, LOW, NUP) + FREQ_FACT * XJCLK_LTE

!           FREE - BOUND:
            ARR(DP, NUP, LOW) = ARR(DP, NUP, LOW) + FREQ_FACT * 
     $                          (PLANCK_FUNC(WAVENUM * light_speed, TL) * (1.0D0 - EXFAC) + XJCLK * EXFAC)

            IF (ARR_LTE_COND) ARR_LTE(DP, NUP, LOW) = ARR_LTE(DP, NUP, LOW) + FREQ_FACT * 
     $                                                (PLANCK_FUNC(WAVENUM * light_speed, TL) *
     $                                                (1.0D0 - EXFAC) + XJCLK_LTE * EXFAC)

          ENDDO L2

!         FINISHING THE CALCULATION OF THE NET RADIATIVE BRACKET
          RRATE(LOW, NUP) = 0.0D0
          RRATE(NUP, LOW) = PI8 * REC * ENLTE(LOW) / ENLTE(NUP)

!         FINISHING THE CALCULATION OF THE ABSOLUTE RADIATIVE RATES
          ARR(DP, LOW, NUP) = C3 * ARR(DP, LOW, NUP)! * 2.0D0 ! The factor of 2 which seemingly should be there due to the trapezoidal integration rule is missing from these formulas
                                                      ! (look into module FGRID for FWEIGHT) but the result seems to be correct so apparently I misunderstood something.

          ARR(DP, NUP, LOW) = C3 * (ENLTE(LOW) / ENLTE(NUP)) * ARR(DP, NUP, LOW)! * 2.0D0

	  IF (ARR_LTE_COND) ARR_LTE(DP, LOW, NUP) = C3 * ARR_LTE(DP, LOW, NUP)

	  IF (ARR_LTE_COND) ARR_LTE(DP, NUP, LOW) = C3 * (POP_LTE(DP, LOW) / POP_LTE(DP, NUP)) * ARR_LTE(DP, NUP, LOW)

          CYCLE LLOW

          !***  LEVELS BELONG TO DIFFERENT ELEMENTS,
          !***  CHARGE DIFFERENCE NOT 1 OR UPPER LEVEL NO GROUND STATE: ZERO RATE

  14      RRATE(LOW, NUP) = 0.0D0
          RRATE(NUP, LOW) = 0.0D0

        ENDDO LLOW

      ENDDO

      !***  ENDLOOP  ---------------------------------------------------------
     
      !***  DIAGONAL ELEMENTS ARE SET TO ZERO

      FORALL (J = 1 : N) RRATE(J, J) =   0.0D0
      FORALL (J = 1 : N) ARR(DP, J, J) = 0.0D0

      RBR(DP, 1 : N, 1 : N) = RRATE(1 : N, 1 : N)

      IF (ARR_LTE_COND) FORALL (J = 1 : N) ARR_LTE(DP, J, J) = 0.0D0

      RETURN

      END SUBROUTINE

      END MODULE
