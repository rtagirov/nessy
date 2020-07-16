      module mod_colli

      use utils

      contains

      subroutine colli(N, ENLTE, TL, ENE, NCHARG, ELEVEL, EINST, CRATE,
     $                 EION, COCO, KEYCOL, WEIGHT, ALTESUM, NATOM, NOM, KODAT,
     $                 levatnum, POPHIIL, POPHML, POPHIL, LEVEL, JOBNUM, DP)

!     collisional transition rates stored in matrix CRATE
!     bound-bound: depending on the element (HE, H)
!     bound-free: not depending on the element (HE, H)
!     called by coma

!     enlte: lte population number for each level

      use hyd_col_rate
      use phys
      use math
      use common_block

      implicit real*8(a - h, o - z)

      integer,intent(in)  :: N, NCHARG, NATOM, NOM, KODAT
      real*8, intent(in)  :: ENLTE,TL,ENE,ELEVEL,EINST,EION,COCO,WEIGHT
      real*8, intent(in)  :: ALTESUM

      integer, intent(in), dimension(N) :: levatnum

      REAL*8, INTENT(IN)  :: POPHIIL, POPHML, POPHIL

      REAL*8, INTENT(OUT) :: CRATE(N, N)

      DIMENSION EINST(N, N)
      DIMENSION ENLTE(N), NCHARG(N), ELEVEL(N), WEIGHT(N)
      DIMENSION EION(N), ALTESUM(4, N)
      DIMENSION COCO(N, N, 4)
      DIMENSION NOM(N)
      DIMENSION KODAT(10)
      CHARACTER*4 KEYCOL(N, N)

!     C1 = H * C / K
      DATA C1 /1.4388D0/ ! (CGS, cm / K)

      REAL*8, DIMENSION(10) :: Temp

      DATA Temp /3D3, 5D3, 7D3, 1D4, 2D4, 3D4, 6.5D4, 1D5, 1.5D5, 2D5/

      INTEGER, INTENT(IN) :: JOBNUM

      INTEGER, INTENT(IN) :: DP

      REAL*8 :: CRJohnUite, CRJohnNorm, CRJeff, RatioJJ, RatioUJ, NewRate

      REAL*8, DIMENSION(6) :: CRUite

      REAL*8 :: ElecProtHydTerm

      LOGICAL :: TRAN_COND, DEPTH_COND

      CHARACTER*10, DIMENSION(N), INTENT(IN) :: LEVEL

      IF (ALLOCATED(ACR)) ACR(DP, 1 : N, 1 : N) = 0.0D0

CMH	collisional electron detachmant of Hminus
CMH	CE: cross section for collisions of Hminus with electrons
CMH	CP: cross section for collisions of Hminus with protons

      THETA = 5040./TL
      THETART = SQRT(THETA)
      THETA32 = THETA*THETART
      THETA13 = THETA**(1./3.)

!      CE = 10*(10.**(-8.7))/THETA32
!      CP = 10*10.**(-7.4)*THETA13

      CE = (10.**(-8.7))/THETA32
      CP = 10.**(-7.4)*THETA13
      CH = 10.**(-10.9+0.5*THETA)/THETA32

      if (natom.gt.30) stop "colli - natom"
      TROOT=SQRT(TL)
      T32=TL*TROOT

      IF (JOBNUM .EQ. 0) THEN

      OPEN(UNIT = 357, FILE = 'Uitenbroek_H_CE.dat', ACTION = 'READ')
      OPEN(UNIT = 358, FILE = 'Uitenbroek_H_CI.dat', ACTION = 'READ')

      OPEN(UNIT = 457, FILE = 'col_rates.out', ACTION = 'WRITE', STATUS = 'REPLACE')

      DO m = 1, 10; WRITE(457, '(13x,F7.0,$)') Temp(m); ENDDO

      WRITE(457, '(/)')

      DO nl = 0, 8

         DO nu = nl + 1, 9

          WRITE(457,'(I2,2x,I2,$)') nu, nl

          x = 1.0D0 - (ISQ(nl + 1) / ISQ(nu + 1))

          DO m = 1, 10

             CRJohnNorm = HYDCOLRATE(nl + 1, nu + 1, Temp(m))
             CRJohnUite = CRJohnNorm * DEXP(x * DABS(HYD_LEV_ENERGY(nl + 1)) / boltz / Temp(m)) / DSQRT(Temp(m))

             WAVENUM = ELEVEL(nu + 2) - ELEVEL(nl + 2)
             WN2 = WAVENUM * WAVENUM
             
             T32 = Temp(m) * DSQRT(Temp(m))

             CRJeff = 3.24D0 * EINST(nu + 2, nl + 2) / WN2 / T32 / (C1 * WAVENUM / Temp(m))**1.68D0

             CRJeff = CRJeff * ENLTE(nu + 2) / ENLTE(nl + 2)

             IF (m .LE. 6) THEN

                IF (m .EQ. 1) READ(357, '(6(ES9.3,2x))') CRUite(1), CRUite(2), CRUite(3), CRUite(4), CRUite(5), CRUite(6)

                CRUite(m) = CRUite(m) * 1D6

                RatioUJ = CRUite(m) / CRJohnUite

             ENDIF

             RatioJJ = CRJeff / CRJohnNorm

             IF (m .LE. 6)                 WRITE(457, '(4x,ES7.1,2x,ES7.1,$)') RatioUJ,  RatioJJ
             IF (m .GT. 6 .AND. m. NE. 10) WRITE(457, '(4x,A7,2x,ES7.1,$)')   '   -   ', RatioJJ
             IF (m .EQ. 10)                WRITE(457, '(4x,A7,2x,ES7.1)')     '   -   ', RatioJJ

          ENDDO

         ENDDO

      ENDDO

      WRITE(457, '(/)')

      nu = 10

      DO nl = 0, 9

          WRITE(457,'(I2,2x,I2,$)') nl, nu

          DO m = 1, 10

             CRJohnNorm = HYDCOLRATE(nl + 1, nu + 1, Temp(m))
             CRJohnUite = CRJohnNorm * DEXP(DABS(HYD_LEV_ENERGY(nl + 1)) / boltz / Temp(m)) / DSQRT(Temp(m))

             TROOT = DSQRT(Temp(m))
      
             G = .1
             EDGE = EION(nl + 2) - ELEVEL(nl + 2)
             EXPFAC = EXP(-C1 * EDGE / Temp(m))
             CRJeff = G * 1.08E-5 * TROOT * EINST(nl + 2, nu + 2) * EXPFAC / EDGE

             IF (ALTESUM(1, nl + 2) .GT. 0.0D0) THEN
                X = 1000.0D0 / Temp(m)
                FOFT = (ALTESUM(3, nl + 2) * X + ALTESUM(2, nl + 2)) * X
                FOFT = 10.0D0**FOFT
                AOFT = ALTESUM(1, nl + 2) * FOFT
                OMSUM = 4.06D0 / TROOT * EXPFAC / EDGE / EDGE / EDGE * AOFT
                CRJeff = CRJeff + OMSUM
             ENDIF

             IF (m .LE. 6) THEN

                IF (m. EQ. 1) READ(358, '(6(ES9.3,2x))') CRUite(1), CRUite(2), CRUite(3), CRUite(4), CRUite(5), CRUite(6)

                CRUite(m) = CRUite(m) * 1D6

                RatioUJ = CRUite(m) / CRJohnUite

             ENDIF

             RatioJJ = CRJeff / CRJohnNorm

             IF (m .LE. 6)                 WRITE(457, '(4x,ES7.1,2x,ES7.1,$)') RatioUJ,  RatioJJ
             IF (m .GT. 6 .AND. m .NE. 10) WRITE(457, '(4x,A7,2x,ES7.1,$)')   '   -   ', RatioJJ
             IF (m .EQ. 10)                WRITE(457, '(4x,A7,2x,ES7.1)')     '   -   ', RatioJJ

          ENDDO

      ENDDO

      CLOSE(357)
      CLOSE(358)
      CLOSE(457)

      ENDIF

      CRATE(1 : N, 1 : N) = 0D0

C***  LOOP OVER ALL TRANSITIONS  ---------------------------------------
      DO 1 NUP = 2, N
	
	DO 1 LOW = 1, NUP - 1

      IF (NOM(LOW) .NE. NOM(NUP)) GOTO 14
      IF (NCHARG(LOW).NE.NCHARG(NUP)) GOTO 8
     
C***  LINE TRANSITION   *******************************************************
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      WN2=WAVENUM*WAVENUM
      WN3=WN2*WAVENUM
C***  ******************************************************************     
C***  HELIUM  ==========================================================

!      IF (NOM(LOW) .EQ. KODAT(1)) THEN
      IF (levatnum(LOW) .EQ. 2) THEN

C***  HE I  *******************************************************************
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCOL  *********
       IF (NCHARG(LOW).EQ.0) THEN
C***  'JEFF': OPTICALLY PERMITTED  TRANSITIONS. JEFFERIES P. 118 (EQ. 6.24)
         IF (KEYCOL(NUP,LOW) .EQ. 'JEFF') THEN
C                                  ====
              OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68
     
C***  'BFK ': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N=2, N=1 OR WITHIN N=2:
C***          BERRINGTON, FON + KINGSTON ('BFK') 1982, MNRAS 200, 347
         ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BFK1') THEN
C                                       ====
              PBFK=COCO(NUP,LOW,1)/TROOT
     +                          +COCO(NUP,LOW,2)+COCO(NUP,LOW,3)*TROOT
     +                          +COCO(NUP,LOW,4)*T32
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
     
         ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BFK2') THEN
C                                       ====
              PBFK=COCO(NUP,LOW,1)/TL
     +                          +COCO(NUP,LOW,2)/TROOT+COCO(NUP,LOW,3)
     +                          +COCO(NUP,LOW,4)*TL
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
     
C***  'BKMS' OR 'BKGR':
C***            FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N=1 OR N.GT.3, N.NE.1:
C***        BENSON + KULANDER ('BK..') 1972, SOLAR PHYSICS 27, 305 (FORMULA 3)
C***        '..MS': MIHALAS + STONE (REF. 11)
C***        '..GR': GREEN (REF. 7)
C***        ATTENTION: COCO(.,.,3) := 1.-ALPHA !!!!!!!!!!!!!!!!!
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BKMS'
C                                        ====
     $        .OR. KEYCOL(NUP,LOW) .EQ. 'BKGR') THEN
C                                        ====
              OMEGA=COCO(NUP,LOW,1)*TL**COCO(NUP,LOW,2)*
     *               EXP(COCO(NUP,LOW,3)*C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
     
C***  'UPS.': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N.NE.1:
C***        EFFECTIVE COLLISIONAL STRENGTH (UPSILON) FROM WERNER SCHMUTZ
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS0') THEN
C                                        ====
            OMEGA=.0
     
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS1') THEN
C                                        ====
              UPSILON=0.05
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
     
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS2') THEN
C                                        ====
              UPSILON=1.
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT

C***  'NONE': TRANSITIONS (FROM N.GT.4) WITH UNKNOWN COLLISIONAL COEFFICIENTS
C***          (COLLISIONAL CROSS SECTION SIGMA(LOW,UP) IS SET TO  PI*A0**2)
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'NONE') THEN
C                                        ====
              OMEGA=5.465E-11*TROOT*(1.+C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
     
C***  'ZERO': NO COLLISIONAL TRANSITION
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'ZERO') THEN
C                                        ====
              OMEGA=0.
     
          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCOL' DECODED

              write (6,*) ' HE I: WRONG KEYWORD KEYCOL '
			STOP 'ERROR'
          ENDIF
C***  HE II  ******************************************************************
       ELSE IF (NCHARG(LOW).EQ.1) THEN
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 119 (EQ. 6.25)
          OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
C***  HE III  *****************************************************************
       ELSE
          write (6,*) '  HE III - LINE TRANSITION  '
          STOP 'ERROR'
       ENDIF
C***  ******************************************************************     
C***  HYDROGEN  ========================================================

!      ELSE IF (NOM(LOW) .EQ. KODAT(2)) THEN
      ELSE IF (levatnum(LOW) .EQ. 1) THEN

C***  H I  *************************************************************
       IF (NCHARG(LOW) .EQ. 0) THEN

C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 118 (EQ. 6.24)
          OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68

!     OMEGA (UP-LOW) according to L. C. Johnson, "Approximations For Collisional
!     And Radiative Transition Rates In Atomic Hydrogen", 1972, ApJ, 174 : 227 - 236
          OMEGA = HYDCOLRATE(LOW - 1, NUP - 1, TL) * ENLTE(LOW) / ENLTE(NUP)

C***	H-MINUS  ************************************************************
       ELSE IF (NCHARG(LOW) .EQ. -1) THEN
        print *,'colli: negative hydrogen found'
        stop
C***  H II  ************************************************************
       ELSE
          write (6,*) 'H II - LINE TRANSITION'
          STOP 'ERROR'
       ENDIF
C***  ******************************************************************
C***  METALS: INDEX GE 3
      ELSE IF (levatnum(LOW) .GE. 3) THEN
C***  ALL IONIZATION STAGES
C***  UPSILON TEMPERATURE INDEPENDENT
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 118 (EQ. 6.24)
C                                        
         IF (KEYCOL(NUP,LOW) .EQ. 'CARB') THEN
C                                  ====
            UPSILON=COCO(NUP,LOW,1)
            IF (UPSILON.LE.0.) THEN
               write (6,*) 'SI - UPSILON .LE. 0'
               STOP 'ERROR'
            ENDIF
            OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
         ELSE
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 119 (EQ. 6.25)
            OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
         ENDIF
C***  *****************************************************
	ELSE

          write(6, *) 'LINE TRANSITION OF UNKNOWN ELEMENT'
          PRINT*,     'LINE TRANSITION OF UNKNOWN ELEMENT'
          PRINT*,     'THE CODE WANTS HELIUM TO BE THE SECOND ELEMENT IN THE DATOM FILE'

          STOP 'ERROR'

      ENDIF

C***  *************************************************************************
C***  COLLISION RATE COEFFICIENTS CRATE
CMH   FOR LINE TRANSITIONS

      CRATE(NUP, LOW) = ENE * OMEGA

!=============================================================================
!     2 - > 3 (LyAlpha) COLLISIONS TEST:

!      TRAN_COND = LOW .EQ. 2 .AND. NUP .EQ. 3

!      DEPTH_COND = (DP .GE. 60) .OR. (DP .GE. 42 .AND. DP .LE. 52)

!      IF (TRAN_COND .AND. DEPTH_COND) CRATE(NUP, LOW) = CRATE(NUP, LOW) * 1D+10
!=============================================================================

      CRATE(LOW, NUP) = CRATE(NUP, LOW) * ENLTE(NUP) / ENLTE(LOW)

      GOTO 1
     
    8 CONTINUE
C***  COLLISIONAL IONIZATION   ************************************************
C***  CHARGE DIFFERENCE MUST BE 1
      IF (NCHARG(NUP) .NE. NCHARG(LOW)+1 ) GOTO 14
C***  UPPER LEVEL MUST BE A GROUND STATE
      IF (NCHARG(NUP) .NE. NCHARG(NUP-1)+1) GOTO 14
C***  OMEGA(LOW-UP) ACCORDING TO JEFFERIES P. 121 (EQ. 6.39)
C***  G = FACTOR DEPENDING ON THE CHARGE OF THE ION (=UPPER LEVEL)
      G=.3
      IF (NCHARG(NUP) .EQ. 2) G=.2
      IF (NCHARG(NUP) .EQ. 1) G=.1
      EDGE = EION(LOW)-ELEVEL(LOW)
      EXPFAC=EXP(-C1*EDGE/TL)
      OMEGA=G*1.08E-5*TROOT*EINST(LOW,NUP)*EXPFAC/EDGE
C***  ADD THE COLLISIONS  TO ADDITIONAL UPPER LEVELS WHICH ARE ASSUMED IN LTE.
C***  NOTE: THE SUMMATION HAS BEEN PERFORMED IN ADVANCE, AND IS IMPLICITELY
C*** CONTAINED IN THE PARAMETERS 'ALTESUM' ENTERING NOW FROM THE ATOMIC DATA
C***  FILE.
      IF (ALTESUM(1,LOW) .GT. .0) THEN

         X=1000./TL
         FOFT=(ALTESUM(3,LOW)*X+ALTESUM(2,LOW))*X
         FOFT=10.**FOFT
         AOFT=ALTESUM(1,LOW)*FOFT
         OMSUM=4.06/TROOT*EXPFAC/EDGE/EDGE/EDGE * AOFT
         OMEGA=OMEGA+OMSUM

      ENDIF

!     OMEGA (LOW-UP) according to L. C. Johnson, "Approximations For Collisional
!     And Radiative Transition Rates In Atomic Hydrogen", 1972, ApJ, 174 : 227 - 236
      IF (LEVEL(NUP) .EQ. 'H II......') OMEGA = HYDCOLRATE(LOW - 1, NUP - 1, TL)

CMH   print a warning if negative values
      IF (ENE .LT. 0.0 .OR. POPHIIL .LT. 0.0 .OR. ENE * OMEGA .LT. 0.0) PRINT*, 'colli: warning!: ', TL, ENE, POPHIIL, cp, ce, OMEGA

C***  BOUND-FREE COLLISIONAL RATE COEFFICIENTS: CRATE
      CRATE(LOW, NUP) = ENE * OMEGA

      CRATE(NUP, LOW) = CRATE(LOW, NUP) * ENLTE(LOW) / ENLTE(NUP)

CMH****************************************************************
CMH	 Transition from negative hydrogen to neutral hydrogen
CMH  LOW = 1: negative hydrogen
CMH	 UP  = 2:  neutral hydrogen
!RT  three processes are important for this transition
CMH	 H- + e- = H + e- + e- (CE * ENE)
CMH	 H- + p+ = H + H       (CP * POPHIIL)
!RT  H- + H  = H + H + e   (CH * POPHIL)
CMH  POPHIIL: population number of protons (HII)
!RT  POPHIL:  population number of neutral hydrogen (HI)
CMH****************************************************************

      Hminfac = 1.0D0

      IF (LEVEL(LOW) .EQ. 'H MINUS..1' .AND. LEVEL(NUP) .EQ. 'H I......1') THEN

         ElecProtHydTerm = Hminfac * (CE * ENE + CP * POPHIIL + CH * POPHIL)
!         ElecProtHydTerm = Hminfac * (CE * ENE + CP * POPHIIL)

         CRATE(LOW, NUP) = CRATE(LOW, NUP) + ElecProtHydTerm
         CRATE(NUP, LOW) = CRATE(NUP, LOW) + ElecProtHydTerm * ENLTE(LOW) / ENLTE(NUP)

CMH****************************************************************
CMH	Transition from protons (H II) to neutral hydrogen
CMH  LOW =  2: neutral hydrogen, ground level
CMH  UP  = 12: protons
CMH	 POPHML: POPNUM of Hminus
CMH	 Hminus collision rates with protons:
CMH	 H- + p+ = H + H (CP * POPHML)
CMH****************************************************************

      ELSEIF (LEVEL(LOW) .EQ. 'H I......1' .AND. LEVEL(NUP) .EQ. 'H II......') THEN

!     The ionization of neutral hydrogen

         HminTerm = Hminfac * CP * POPHML

         CRATE(NUP, LOW) = CRATE(NUP, LOW) + HminTerm

CMH	transition up = 12 to low=2 is
CMH	1. the recombination of protons to neutral hydrogen and
CMH	2. the ionization of negative hydrogen due to proton collision

         CRATE(LOW, NUP) = CRATE(LOW, NUP) + HminTerm * ENLTE(NUP) / ENLTE(LOW)

      ENDIF

      GOTO 1

C***  LEVELS BELONG TO DIFFERENT ELEMENTS,
C***  CHARGE DIFFERENCE NOT 1  OR UPPER LEVEL NO GROUND STATE: ZERO RATE

   14 CRATE(NUP, LOW) = 0D0; CRATE(LOW, NUP) = 0D0

    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
     
C***  DIAGONAL ELEMENTS ARE SET TO ZERO
      FORALL (J = 1 : N) CRATE(J, J) = 0D0

      IF (ALLOCATED(ACR)) ACR(DP, 1 : N, 1 : N) = CRATE(1 : N, 1 : N)

      RETURN

      END SUBROUTINE

      END MODULE
