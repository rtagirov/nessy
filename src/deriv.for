      MODULE MOD_DERIV

      CONTAINS

      SUBROUTINE DERIV(DLOWUP,I,NUP,LOW,IND,
     $                 NPLUS1,EN,CRATE,RRATE,EXPFAC,NFEDGE,
     $                 WCHARM,ND,L,TL,ENLTE,PHI,PWEIGHT,NFL,DELTAX,XMAX,
     $                 DETAL,DOPAL,SLNEW,OPAL,XJLAPP,XJCAPP,
     $                 FWEIGHT,DOPA,DETA,OPAC,SCNEW,XLAMBDA,NF,
     $                 N,NCHARG,WEIGHT,ELEVEL,NOM,EINST,SIGMAKI,LASTIND,
     $                 XJL, AccFact, SLOLD)

C*******************************************************************************
C***  DERIVATIVE OF RATE EQ. COEFFICIENTS :
C***  D(I,M,J):= N(M) * D(RATCO(M,J))/DN(I) - N(J) * D(RATCO(J,M))/DN(I)
C*******************************************************************************

      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      REAL*8, INTENT(IN), DIMENSION(N, N) :: CRATE ! Collision Rate

      REAL*8, dimension(N)      :: ENLTE  ! population number in LTE
      REAL*8, dimension(NPLUS1) :: EN     ! population number in NLTE
      REAL*8, dimension(ND, NF) :: WCHARM

      DIMENSION NCHARG(N),RRATE(N,N)
      DIMENSION EINST(N,N),ELEVEL(N)
      DIMENSION WEIGHT(N),NFEDGE(N)
      DIMENSION NOM(N)
      DIMENSION DOPA(NF),DETA(NF),OPAC(NF),SCNEW(NF)
      DIMENSION FWEIGHT(NF),XLAMBDA(NF),XJCAPP(NF),EXPFAC(NF)
      DIMENSION DOPAL(LASTIND),SLNEW(LASTIND),OPAL(LASTIND),DETAL(LASTIND)
      DIMENSION XJLAPP(LASTIND)
      DIMENSION PHI(NFL),PWEIGHT(NFL)
      DIMENSION SIGMAKI(NF,N)

      REAL*8, DIMENSION(LASTIND), INTENT(IN) :: XJL, AccFact, SLOLD

C***  C2 = 2 * H * C    ( H AND C IN CGS UNITS )
      DATA C2 / 3.9724d-16 /
C***  C3 = 4 * PI / H / C
      DATA C3 / 0.06327d18 /
C***  PI8 = 8 * PI

      DATA PI8 / 25.13274123d0 /

      DLOWUP = 0.0D0

      IF (NOM(NUP) .NE. NOM(LOW)) RETURN

      NDCH = NCHARG(NUP) - NCHARG(LOW)

      IF (NDCH .GT. 1) RETURN

      IF (NDCH .EQ. 0) GOTO 1

C***  M - J IS BOUND-FREE TRANSITION  ******************************************
C***  UPPER LEVEL MUST BE A GROUND STATE

      IF (NCHARG(NUP) .EQ. NCHARG(NUP - 1)) RETURN

C***  DERIVATIVE WITH RESPECT TO ELECTRON DENSITY
      IF (I .EQ. NPLUS1) THEN

         DLOWUP = - (EN(LOW) * CRATE(LOW, NUP) - EN(NUP) * (2.d0 * CRATE(NUP, LOW) + RRATE(NUP, LOW))) / EN(NPLUS1)

         IF (CONST_ELEC) DLOWUP = 0.0D0 ! PRE-SET ELECTRON CONCENTRATION

      ENDIF

      if (isnan(dlowup)) stop 'deriv bf elec is NaN'

C***  RATE INTEGRAL (DERIVATIVE OF NET RADIATIVE BRACKET)

      SUM = 0.0D0

C***  PREFACTOR FOR THE DERIVATIVE OF 1/S-BOUND-FREE WITH RESPECT TO N(I)

      IF (I .EQ. LOW) THEN

         PREFAC = -1.0D0 / C2 / EN(NUP)

      ELSEIF (I .EQ. NUP) THEN

         PREFAC = EN(LOW) / EN(NUP) / EN(NUP) / C2

      ELSEIF (I .EQ. NPLUS1) THEN

         PREFAC = EN(LOW) / EN(NUP) / EN(NPLUS1) / C2

      ELSE

         PREFAC = 0.0D0

      ENDIF

      NFLOW = NFEDGE(LOW)

      DO K = 1, NFLOW

         WAVENUM = 1.0D8 / XLAMBDA(K)

         W2 = WAVENUM * WAVENUM

         W3 = W2 * WAVENUM

C***  CALCULATE BOUND-FREE SOURCE FUNCTION FOR TRANSITION LOW-UP ONLY

         EXFAC = EXPFAC(K)

         G = EXFAC * ENLTE(LOW) / ENLTE(NUP)

         SBF = C2 * W3 / (EN(LOW) / (EN(NUP) * G) - 1.0D0)

         DSUM = XJCAPP(K) * PREFAC / W3 / G + WCHARM(L, K) * (DOPA(K) * SCNEW(K) - DETA(K)) / (OPAC(K) * SBF)

C***  MAKE USE OF SIGMAKI, THE TABULATED PHOTOCROSSECTIONS IN CM**2

         SUM = SUM + DSUM * SIGMAKI(K, LOW) * EXFAC * W2 * FWEIGHT(K)

      ENDDO

      DLOWUP = DLOWUP + SUM * PI8 * EN(NUP) * ENLTE(LOW) / ENLTE(NUP)

      if (isnan(dlowup)) then 

         print*, 'deriv bf dlowup is NaN'

         stop

      endif

      RETURN

C***  M - J IS BOUND-BOUND TRANSITION   ****************************************
    1 CONTINUE

C***  DERIVATIVE WITH RESPECT TO ELECTRON DENSITY
      IF (I .EQ. NPLUS1) THEN

          CRRatio = EN(LOW) * CRATE(LOW, NUP) / EN(NUP) / CRATE(NUP, LOW)

          DLOWUP = EN(NUP) * CRATE(NUP, LOW) * (1.0D0 - CRRatio) / EN(NPLUS1)
!         DLOWUP = (EN(NUP) * CRATE(NUP, LOW) - EN(LOW) * CRATE(LOW, NUP)) / EN(NPLUS1)

         IF (CONST_ELEC) DLOWUP = 0.0D0 ! PRE-SET ELECTRON CONCENTRATION

          RETURN

      ENDIF

      if (isnan(dlowup)) stop 'deriv bb elec dlowup is NaN'

!     DERIVATIVES WITH RESPECT TO POPNUMBER EN(I)

!     RINAT TAGIROV:
!     These two processes are not taken into account in the Jacobian calculations:
!     H- + e- <-> H + e- + e-
!     H- + p+ <-> H + H
!     Since the convergence was attained we decided not to bother with including them in this procedure.

!     ESTIMATE OF CURRENT LINE INDEX IND
      IND = IND + 1

!     THE DERIVATIVES OF THE CONT. BACKGROUND-OPACITY AND EMISSIVITY
!     ARE NEGLECTED. HENCE:
      IF (I .NE. NUP .AND. I .NE. LOW) RETURN

!     RUDIMENTAL LINES: ZERO CORE ASSUMED
      IF (EINST(LOW, NUP) .EQ. -2.d0) RETURN

!     TERM FOR NET RADIATIVE BRACKETS

      IF (SLNEW(IND) .NE. 0.0D0) THEN

         ETAL = SLNEW(IND) * OPAL(IND)

         DLOWUP_RB = XJL(IND) * (OPAL(IND) * DETAL(IND) - DOPAL(IND) * ETAL) / ETAL**2.0D0

         DLOWUP_RB = EINST(NUP, LOW) * DLOWUP_RB

!     RINAT TAGIROV:
!     Taking the local operator acceleration into account in the Jacobian calculations
!-------------------------------------------------------------------------------------
         DLOWUP_RB = (1.0D0 - AccFact(IND) * SLOLD(IND) / XJL(IND)) * DLOWUP_RB
!-------------------------------------------------------------------------------------

      ELSE

         DLOWUP_RB = 0.0D0

      ENDIF

      DLOWUP = EN(NUP) * DLOWUP_RB

      if (isnan(dlowup)) stop 'deriv bb dlowup is NaN'

      RETURN

      END SUBROUTINE

      END MODULE
