      module mod_linpop

      contains

      SUBROUTINE LINPOP(T,RNE,ENTOT,ITNE,POPNUM,DEPART_ZWAAN,POP1,
     $                  N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EINST,LEVEL,
     $                  XLAMBDA,FWEIGHT,XJC,NF,XJL,WCHARM,EPSILON,NODM,
     $                  DELTAC,MODHEAD,JOBNUM,IFRRA,ITORA,
     $                  RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $                  VELO,GRADI,VDOP,INDNUP,INDLOW,LASTIND,
     $                  OPAC,DOPA,DETA,SIGMAKI,
     $                  ND,LSRAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,
     $                  ALTESUM,ETAC,NFEDGE,EXPFAC,NOM,NATOM,KODAT,NFIRST,
     $                  NLAST,WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,JOBMAX)

!     see Koesterke et al 1992 A&A 255, 490

!     CALCULATION OF NEW NLTE POPULATION NUMBERS (ARRAY POPNUM)
!     RADIATIVE RATES ARE CALCULATED WITH THE SCHARMER RADIATION FIELD
!     THE (HENCE NON-LINEAR) RATE EQUATIONS ARE SOLVED BY LINEARIZATION
!     POP1 = OLD POPULATION NUMBERS
!     POPNUM = NEW POPULATION NUMBERS
!     RNE = REL. ELECTRON DENSITY - UPDATED IN THIS SUBROUTINE
!     EN(J) = NEW NLTE POP. NUMBERS  AT CURRENT DEPTH POINT
!     ENLTE(J) = LTE POP. NUMBERS AT CURRENT DEPTH POINT
!     ENDELTA(J) = DELTA OF NEW - CURRENT POP. NUMBERS
!     DB(N, N) = BROYDEN MATRIX AT CURRENT DEPTH POINT
!     DM(N, N) = ONLY NEEDED IF NO BROYDEN MATRIX FOUND AT START OF NEW MODEL
!     IF THERE IS A DB MATRIX, DM IS USED AS HELP-MATRIX FOR ALGORITHM.
!     VOLD(J) = VEKTOR OF K-1 ITERATION USED FOR K ITERATION
!     V1(J) = B-VEKTOR COMPUTED IN COMA()
!     V2 - V5(J) = HELP-VEKTOR FOR BROYDEN ALGORITHMUS

      USE MOD_BFCROSS
      use MOD_ERF_INF
      USE MOD_ERROR
      use MOD_FLGRID
      use MOD_ISRCHFGT
      use MOD_LTEPOP
      use MOD_PRIRAT
      USE MOD_XRUDI
      use MOD_CCORE
      use MOD_COMA
      use ABUNDANCES
      use MOD_LIOP

      use file_operations
      use common_block
      use matoper
      use broyden

      IMPLICIT REAL*8(A - H, O - Z)

      integer, intent(in) :: lastind, N

      real*8, dimension(lastind) :: opal

      real*8, dimension(N + 1)        :: V1, V2, ENDELTA, V4, V5, VOLD
      real*8, dimension(N + 1, N + 1) :: DB, DM

      real*8, intent(in),  dimension(N, N) :: einst
      real*8, intent(in),  dimension(N)    :: ELEVEL, EION

      integer, intent(in),  dimension(N)   :: NCHARG, NOM

      integer, intent(out), dimension(N) :: NFEDGE

      real*8, dimension(N, N) :: CRATE, RRATE

      real*8, dimension(N + 1) :: en

      character*8, intent(in), dimension(N) :: agaunt

      DIMENSION T(ND),RNE(ND),ENTOT(ND),ITNE(ND)
      DIMENSION XLAMBDA(NF),EXPFAC(NF)

      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION SCOLD(NF, ND)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      REAL*8 SIGMAKI(NF), ALPHA(N), SEXPO(N)

      real*8, allocatable, dimension(:, :) :: RATCO

      real*8 rnel, rnellte

      real*8, allocatable :: ABXYZ_new(:)

      LOGICAL NODM, NEWRAP, STRONG_CONV, WEAK_CONV

      real*8 :: part1, part2, epsdn

      COMMON / COMNEGI / NEGINTL
      CHARACTER modhead*104
      CHARACTER LEVEL(N)*10
      CHARACTER*4 KEYCOL(N,N)
      DIMENSION WAVARR(N,NF),SIGARR(N,NF)

      REAL*8 :: POPHIIL, POPHML, POPHIL

      real*8, allocatable, dimension(:) :: phi, pweight

      real*8  WEIGHT(*)
      logical LBKG
      integer XLBKG1,  XLBKG2, iii

      character*10, dimension(*) :: MAINPRO,MAINLEV
      integer,      dimension(*) :: IWARN

      integer,      dimension(NATOM) :: KODAT

      real*8,       dimension(ND, NF) :: WCHARM
      real*8,       dimension(NF) ::     FWEIGHT

      real*8,       dimension(*) :: DETA, DOPA

      real*8,       dimension(LASTIND) :: DETAL, DOPAL

      real*8,       dimension(*) :: ETA, ETAC, GRADI
      real*8,       dimension(*) :: OPA, OPAC
      real*8,       dimension(*) :: THOMSON

      real*8,       dimension(N, N, 4) :: COCO
      real*8,       dimension(4, N) ::    ALTESUM
      real*8,       dimension(ND) ::      RADIUS, VELO

      REAL*8, DIMENSION(ND, LASTIND) :: XJL, JNEW
      REAL*8, DIMENSION(ND, NF) ::      XJC

      REAL*8, DIMENSION(LASTIND) :: XJLAPP

      REAL*8, DIMENSION(:), ALLOCATABLE :: XJC_EDG

      REAL*8, DIMENSION(N) ::     ENLTE
      REAL*8, DIMENSION(ND, N) :: POPNUM, POPNUM_LTE, POP1
      REAL*8, DIMENSION(ND, N) :: DEPART_ZWAAN

      REAL*8, DIMENSION(ND, LASTIND) :: SLOLD, SLNEW

      REAL*8 :: opalind, etalind

      REAL*8, DIMENSION(ND) :: ELEC_CONC_OLD, ElecConc, ElecConcLTE, ElecConcDep

      REAL*8, DIMENSION(ND) :: ONE, Z

      INTEGER, INTENT(IN) :: JOBNUM

      INTEGER :: IND, LOW, NUP

      CHARACTER(:), ALLOCATABLE :: CONV_FILE

      LOGICAL :: NRRM_FILE_EXISTS, NCRM_FILE_EXISTS, NTRM_FILE_EXISTS, CONV_FILE_EXISTS

      LOGICAL :: PRINT_LTE_ARR

      REAL*8 :: ONE_PRO, PRO_LEV, TWO_PRO, LEV_TWO, TWO_ONE

      REAL*8, DIMENSION(ND) :: ELEC_CONC, HEAVY_ELEM_CONC

      INTEGER :: N_HI_LEV, N_HI_LIN

C***  C1 = H * C / K (CM * KELVIN)
      DATA C1 /1.4388D0/

      CONV_FILE = CONV_DIR//'ALL'

      PRINT_LTE_ARR = .FALSE.

      DEPART_ZWAAN(1 : ND, 1 : N) = 0.0D0

      JNEW(1 : ND, 1 : LASTIND) = 0.0D0

      ElecConc(1 : ND) = 0.0D0
      ElecConcLTE(1 : ND) = 0.0D0

      ELEC_CONC_OLD(1 : ND) = 0.0D0

      ElecConcDep(1 : ND) = 0.0D0

      ALLOCATE(ARR(ND, N, N)); ARR(:, :, :) = 0.0D0
      ALLOCATE(ACR(ND, N, N)); ACR(:, :, :) = 0.0D0
      ALLOCATE(RBR(ND, N, N)); RBR(:, :, :) = 0.0D0

      Z(1 : ND) =   0.0D0
      ONE(1 : ND) = 1.0D0

      NRANK = N + 1

      NPLUS1 = N + 1

      if (allocated(ratco)) deallocate(ratco); allocate(ratco(nrank, nrank))

      SLOLD(1 : ND, 1 : LASTIND) = 0.0D0

C***  REMOVE NEGATIVE LINE INTENSITIES

      NEGINTL = 0

      DO L = 1, ND

         DO IND = 1, LASTIND

            LOW = INDLOW(IND)

            NUP = INDNUP(IND)

            IF (EINST(LOW, NUP) .NE. -2.0D0) THEN

               IF (XJL(L, IND) .LT. 0.0D0) THEN

                  XJL(L, IND) = ABS(XJL(L, IND))
 
                  NEGINTL = NEGINTL + 1

               ENDIF

               XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

!              CALCULATE THE OPACITY OPAL FOR A GIVEN LINE AND DEPTH POINT
               CALL LIOP_SBE(EINST(NUP, LOW), WEIGHT(LOW), WEIGHT(NUP), LOW, NUP,
     $                       XLAM, ENTOT(L), POP1(L, 1 : N), RSTAR, opalind, etalind, VDOP, N)

               IF (opalind .LE. 0.0D0) THEN

                  SLOLD(L, IND) = 0.0D0

               ELSE

                  SLOLD(L, IND) = etalind / opalind

               ENDIF

            ENDIF

         ENDDO

      ENDDO

!     GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
!     SIGMAKI(K,LOW) IN CM**2
      CALL BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,
     $             XLAMBDA(1 : NF),ALPHA,SEXPO,AGAUNT,NOM,
     $             WAVARR(:, 1 : NF),SIGARR(:, 1 : NF))

!     DETERMINE SCOLD AT ALL DEPTH POINTS
      CALL CCORE(NF,DELTAC,MODHEAD,JOBNUM,
     $           SCOLD,RADIUS,XLAMBDA,ND,T,RNE,POP1,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $           N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,SIGMAKI,
     $           WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2)

!     GENERATE LINE FREQUENCY GRID AND INTEGRATION WEIGHTS
      CALL FLGRID(NFL, PHI, PWEIGHT, DELTAX)

      XMAX = DELTAX * FLOAT(NFL - 1) * .5
      ERXMIN = ERF_INF(-XMAX)

      N_HI_LEV = NLAST(1) - 2

      N_HI_LIN = 0

      DO J = 1, N_HI_LEV - 1; N_HI_LIN = N_HI_LIN + J; ENDDO

      PRINT*, 'NUMBER OF HI LINES: ', N_HI_LIN

      INQUIRE(FILE = CONV_FILE, EXIST = CONV_FILE_EXISTS)

      IF (.NOT. CONV_FILE_EXISTS) LAMBDA_ITER = 0

      IF (CONV_FILE_EXISTS) LAMBDA_ITER = NUM_OF_LINES(CONV_FILE) - 1

      IF (.NOT. LTE_RUN .AND. LAMBDA_ITER .EQ. 0) THEN

         CALL MKDIR(NLTE_DIR_1); CALL CLEAN_DIR(NLTE_DIR_1)
         CALL MKDIR(NLTE_DIR_2); CALL CLEAN_DIR(NLTE_DIR_2)

         DO J = 1, N; CALL PRINT_NLTE_LEV(LEVEL(J), POP1(1 : ND, J), POP1(1 : ND, J), ONE); ENDDO

         CALL PRINT_NLTE_LEV('ELECTRONS ', RNE(1 : ND), RNE(1 : ND), ONE)

         CALL PRINT_HYD_NLTE_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, T,
     $                           Z, Z, Z, Z, Z,
     $                           POP1(1 : ND, 1), POP1(1 : ND, 2),
     $                           Z, Z, Z, Z, Z)

         DO IND = 1, N_HI_LIN
 
            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            CALL PRINT_HYD_NLTE_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, T,
     $                              LLO(1 : ND, IND), 
     $                              XJL(1 : ND, IND),
     $                              SLOLD(1 : ND, IND),
     $                              Z,
     $                              SLOLD(1 : ND, IND),
     $                              POP1(1 : ND, LOW), POP1(1 : ND, NUP),
     $                              Z, Z, Z, Z, Z)

         ENDDO

         DO I = 2, NLAST(1) - 1

            CALL PRINT_HYD_NLTE_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, T,
     $                              Z, Z, Z, Z, Z,
     $                              POP1(1 : ND, I), POP1(1 : ND, NLAST(1)),
     $                              Z, Z, Z, Z, Z)

         ENDDO

         LAMBDA_ITER = 1

         IF (JOBMAX .LE. 5.) STOP 'LTE IS DONE.'

      ELSEIF (LTE_RUN) THEN

         CALL MKDIR(LTE_DIR_1); CALL CLEAN_DIR(LTE_DIR_1)
         CALL MKDIR(LTE_DIR_2); CALL CLEAN_DIR(LTE_DIR_2)

         ALLOCATE(XJL_LTE(ND, LASTIND)); XJL_LTE(1 : ND, 1 : LASTIND) = XJL(1 : ND,  1 : LASTIND)
         ALLOCATE(XJC_LTE(ND, NF));      XJC_LTE(1 : ND, 1 : NF) =      XJC(1 : ND,  1 : NF)
         ALLOCATE(POP_LTE(ND, N));       POP_LTE(1 : ND, 1 : N) =       POP1(1 : ND, 1 : N)

         CALL PRINT_LTE_LEV('ELECTRONS ', 115, RNE(1 : ND))

         DO J = NFIRST(1), NLAST(1); CALL PRINT_LTE_LEV(LEVEL(J), J, POP_LTE(1 : ND, J)); ENDDO

         ALLOCATE(ARR_LTE(ND, N, N)); ARR_LTE(1 : ND, 1 : N, 1 : N) = 0.0D0
         ALLOCATE(XJC_EDG(ND));       XJC_EDG(1 : ND) = 0.0D0

         PRINT_LTE_ARR = .TRUE.

      ENDIF

      ALLOCATE(NOFILE(ND))
 
!     OPEN FORT19 FOR READING THE BROYDEN MATRIX
      CALL DBOPEN(MODHEAD)

      NOFILE(1 : ND) = .TRUE.

!     LOOP OVER ALL DEPTH POINTS
      DO 100 L = 1, ND

      ELEC_CONC_OLD(L) = SUM(NCHARG(1 : N) * POP1(L, 1 : N))

      if (.NOT. allocated(ABXYZ_new)) allocate(ABXYZ_new(NATOM))
      ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)

!     PRE - CALCULATE FREQUENCY INDICES OF IONIZATION EDGES
      DO NA = 1, NATOM

         DO LOW = NFIRST(NA), NLAST(NA) - 1

            EDGE = EION(LOW) - ELEVEL(LOW)

            EDGELAM = 1.0E+8 / EDGE

            NFEDGE(LOW) = ISRCHFGT(NF, XLAMBDA, 1, EDGELAM) - 1

         ENDDO

      ENDDO

!     TRY FIRST BROYDEN ITERATION
      NEWRAP = .FALSE.

!     MAXIMUM NUMBER OF BROYDEN ITERATIONS
      ITMAX = 30
!      ITMAX = 1000

!     DEMANDED ACCURACY OF THE BROYDEN (OR NEWTON) ITERATION
      EPSDN = EPSILON * 2.0D-4
!      EPSDN = EPSILON * 2.0D-8

      IF (EPSDN .GT. 1.0D-5) EPSDN = 1.0D-5
      IF (EPSDN .LT. 1.0D-8) EPSDN = 1.0D-8

!     BROYDEN ITERATION
!     THE ITERATION STARTS FROM THE OLD POPULATION NUMBERS

!     ENTRY POINT FOR FAILED BROYDEN ITERATION
      IRESTA = 0

   50 CONTINUE

      EN(1 : N) = POP1(L, 1 : N)

      EN(NPLUS1) = RNE(L)

      TL = T(L)

!     LOAD BROYDEN-MATRIX OF CURRENT DEPTH POINT
      IF (.NOT. NOFILE(L)) CALL DBLOAD(DB, L, NRANK)

      !BEGINNING OF THE BROYDEN ITERATION LOOP
      ITNE(L) = 0
   10 ITNE(L) = ITNE(L) + 1

!     PRE-CALCULATE EXPONENTIAL FACTORS FOR THE TEMPERATURE OF THE CURRENT DEPTH POINT
!     THIS MUST BE REPEATED, WHEN THE TEMPERATURE HAS BEEN UPDATED
      IF (ITNE(L) .EQ. 1) THEN

          DO K = 1, NF

             WAVENUM = 1.E8 / XLAMBDA(K)

             EXPFAC(K) = EXP(-C1 * WAVENUM / TL)

          enddo

      ENDIF

      !***  CALCULATE LTE POP. NUMBERS
      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ_new,NFIRST,NLAST,NATOM)

!     CALCULATE NODM, NODM = TRUE => DM IS NOT NEEDED
      IF ((ITNE(L) .EQ. 1 .AND. NOFILE(L)) .OR. NEWRAP) THEN; NODM = .FALSE.; ELSE; NODM = .TRUE.; ENDIF

!     SETUP COEFFICIENT MATRICES
      POPHIIL = POPNUM(L, NLAST(1)) * ENTOT(L)
      POPHML  = POPNUM(L, 1) *        ENTOT(L)
      POPHIL =  POPNUM(L, 2) *        ENTOT(L)

!     LLO is declared in comblock.for
      CALL COMA(CRATE,RRATE,RATCO,DM,N,NRANK,V1,ABXYZ_new,
     $          ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,EION,WEIGHT,ALTESUM,
     $          XLAMBDA,FWEIGHT,XJC,NF, L, XJL(L, 1 : LASTIND), ND, XJLAPP,
     $          SLOLD(L, 1 : LASTIND), LASTIND, INDLOW,
     $          INDNUP,NOM,NATOM,KODAT,NFIRST,NLAST,PHI,PWEIGHT,DELTAX,XMAX,
     $          NFL,OPAC,DOPA,DETA,OPAL,SLNEW(L, 1 : LASTIND),
     $          DOPAL, DETAL, SIGMAKI,ETAC,NFEDGE,EXPFAC,NODM,
     $          WCHARM,EN,RSTAR,SCOLD,VDOP,COCO,KEYCOL,
     $          POPHIIL,POPHML, POPHIL, LLO(L, 1 : LASTIND), ITNE(L),
     $          LEVEL, JOBNUM, IRESTA)

      JNEW(L, 1 : LASTIND) = XJLAPP(1 : LASTIND)

      IF (NEWRAP) THEN

!       ALGEBRA OF ONE NEWTON ITERATION STEP
         CALL VMF(V2, EN, RATCO, NRANK)

         CALL VSUB(V1, V2, NRANK) ! V1 = V1 - V2

         CALL INV(NRANK, DM)

         CALL VMF(ENDELTA, V1, DM, NRANK) ! ENDELTA = DM * V1a

      ELSE

!        ALGEBRA OF ONE BROYDEN ITERATION STEP
!        calculate the resulting vector V4 using the current populations EN

         CALL VMF(V4, EN, RATCO, NRANK)

        !*** calculate the error vector V1 = V1 - V4

         CALL VSUB (V1, V4, NRANK)

         IF (ITNE(L) .EQ. 1 .and. nofile(L)) THEN

!           IF FORT.19 NOT FOUND AND THE FIRST ITERATION STEP, DB IS DM^T
            CALL INV(NRANK, DM)

            CALL ACOPY(DB, DM, NRANK)

         ELSE

!       ALGEBRA FOR THE BROYDEN FORMULA
!       delta_y = VOLD = VOLD - V1
            CALL VSUB(VOLD, V1, NRANK)

!       V2 = delta_y * B
            CALL VMF(V2, VOLD, DB, NRANK)

!       V3 = V2 * ENDELTA; ENDELTA is the difference between pops = delta_x; V3 - the denominator (a number)

            CALL VMV(V3, V2, ENDELTA, NRANK)

            V5(1 : NRANK) = ENDELTA(1 : NRANK)

!            do i = 1, nrank; print*, 'endelta after:', i, endelta(i); enddo; stop

!           ENDELTA = B * delta_x = DB * V5
            CALL VMT(ENDELTA, DB, V5, NRANK)

!           V5 = delta_x - delta_y * B = V5 - V2; V5 = V5 - V2
            CALL VSUB(V5, V2, NRANK) 

!           DM = B * delta_x * V5; numerator (matrix) DM = dyadic product ENDELTA * V5
            CALL VMD(DM, ENDELTA, V5, NRANK)

!           matrix DM (numerator) divided by V3
            CALL VDIFF(DM, V3, NRANK)

!           new Broyden matrix DB = DB + DM
            CALL VADDM(DB, DM, NRANK)

         ENDIF

         VOLD(1 : NRANK) = V1(1 : NRANK)

!        calculate the new population-correction vector ENDELTA = V1 * DB
         CALL VMF(ENDELTA, V1, DB, NRANK)

      ENDIF

!     ALGEBRA COMMON TO BOTH APPROACHES; UPDATE THE POPULATIONS EN = EN + ENDELTA
      CALL VADD(EN, ENDELTA, NRANK)

!     CONVERGENCE CHECK FOR BROYDEN OR NEWTON ITERATION
      STRONG_CONV = .TRUE.

      DO J = 1, NPLUS1; STRONG_CONV = STRONG_CONV .AND. (ABS(ENDELTA(J) / EN(J)) .LT. EPSDN .OR. ABS(EN(J)) .LT. 1.0D-15); ENDDO

      IF (.NOT. STRONG_CONV .AND. ITNE(L) .LT. ITMAX) THEN

         GOTO 10

      ELSEIF (.NOT. STRONG_CONV .AND. ITNE(L) .GE. ITMAX) THEN

         WEAK_CONV = .TRUE.

         DO J = 1, NPLUS1

            WEAK_CONV = WEAK_CONV .AND. (ABS(ENDELTA(J) / EN(J)) .LT. 1.0D+2 * EPSDN .OR. ABS(EN(J)) .LT. 1.0D-10)

         ENDDO

         IF (WEAK_CONV) THEN

            WRITE(*, '(/,90x,A)'), 'ONLY WEAK CONVERGENCE HAS BEEN REACHED'

         ELSE

            WRITE(*, '(/,A,/)'), 'CONVERGENCE HAS NOT BEEN REACHED...'

            IRESTA = IRESTA + 1

            IF (IRESTA .GE. 2) THEN

               WRITE(*, '(/,A,/)'), 'NO SOLUTION FOUND...'

               ITNE(L) = -ITMAX

            ELSE

               WRITE(*, '(A)'), 'NEWTON - RAPHSON SCHEME WILL BE IMPLEMENTED...'

               NOFILE(L) = .TRUE.

               NEWRAP = .TRUE.

               ITMAX = 9

               EPSDN = EPSILON * 1D-1

               GOTO 50

            ENDIF

         ENDIF

      ENDIF

!     END OF BROYDEN/NEWTON-ITERATION LOOP

!     SAVE BROYDEN-MATRIX IN FORT19
      IF (NEWRAP) CALL ACOPY(DB, DM, NRANK)

      CALL DBSAVE(DB, L, NRANK)

!     updating array popnum and calculating the departure coefficients:
!     Rinat Tagirov:
!     the most accurate way to do this is to update the lte populations here because upon the last
!     call of the ltepop procedure the electron concentration has not
!     yet converged, i.e. if we don't do that here the lte populations will
!     correspond to the penultimate broyden/newton iteration
!*******************************************************************************************

      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ_new,NFIRST,NLAST,NATOM)
!*******************************************************************************************

      POPNUM(L, 1 : N) =       EN(1 : N)
      POPNUM_LTE(L, 1 : N) =   ENLTE(1 : N)
      DEPART_ZWAAN(L, 1 : N) = EN(1 : N) / ENLTE(1 : N)

C***  UPDATING THE ELECTRON DENSITY
      RNE(L) = EN(NPLUS1)

C***  UPDATING THE ELECTRON TEMPERATURE

      ElecConc(L) =    SUM(NCHARG(1 : N) * EN(1 : N))
      ElecConcLTE(L) = SUM(NCHARG(1 : N) * ENLTE(1 : N))

      ElecConcDep(L) = ElecConc(L) / ElecConcLTE(L)

C***  PRINTOUT OF RATE COEFFICIENTS ETC.  ------------------------------
      IF (LSRAT.NE.-1) THEN
      IF ((L.GE.IFRRA.AND.L.LE.ITORA).OR.ITORA.EQ.0) THEN
      NETTO=1
      LM1=L-1
      IF (IFRRA.GT.0) LM1=L-IFRRA
      IF  (((LM1)/LSRAT)*LSRAT.EQ.(LM1).OR.L.EQ.ND)
     $CALL       PRIRAT (ITNE(L),N,LEVEL,L,CRATE,RRATE,RATCO,EN,
     $           IFRRA,MODHEAD,JOBNUM,NETTO )
      ENDIF
      ENDIF
      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)
  100 CONTINUE

!     CLOSE FORT19 FOR BROYDEN-MATRIX
      CALL DBCLOSE(ND, NRANK, DB)

!     CALCULATING CORMAX FOR ELECTRONS
      CORMAX_ELEC = MAXVAL(DABS(1.0D0 - ELEC_CONC_OLD(3 : ND - 1) / ElecConc(3 : ND - 1)))

      IF (.NOT. LTE_RUN) THEN

!         CALL PRINT_NLTE_LEV('H MINUS..1', NFIRST(1), POPNUM_LTE(1 : ND, NFIRST(1)),
!     $                       POPNUM(1 : ND, NFIRST(1)), DEPART_ZWAAN(1 : ND, NFIRST(1)), ONE)

         DO J = 1, N; CALL PRINT_NLTE_LEV(LEVEL(J), POPNUM_LTE(1 : ND, J),
     $                                    POPNUM(1 : ND, J), DEPART_ZWAAN(1 : ND, J)); ENDDO

!         CALL PRINT_NLTE_LEV('H II......', NLAST(1), POPNUM_LTE(1 : ND, NLAST(1)),
!     $                       POPNUM(1 : ND, NLAST(1)), DEPART_ZWAAN(1 : ND, NLAST(1)), ONE)

         CALL PRINT_NLTE_LEV('ELECTRONS ', ElecConcLTE(1 : ND),
     $                        ElecConc(1 : ND), ElecConcDep(1 : ND))

         CALL PRINT_HYD_NLTE_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, T,
     $                           Z, Z, Z, Z, Z,
     $                           POPNUM(1 : ND, NFIRST(1)), 
     $                           POPNUM(1 : ND, NFIRST(1) + 1),
     $                           ARR(1 : ND, NFIRST(1) + 1, NFIRST(1)), ARR(1 : ND, NFIRST(1), NFIRST(1) + 1),
     $                           ACR(1 : ND, NFIRST(1) + 1, NFIRST(1)), ACR(1 : ND, NFIRST(1), NFIRST(1) + 1),
     $                           RBR(1 : ND, NFIRST(1) + 1, NFIRST(1)))

         DO IND = 1, N_HI_LIN
 
            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            CALL PRINT_HYD_NLTE_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, T,
     $                              LLO(1 : ND, IND), 
     $                              XJL(1 : ND, IND), 
     $                              SLOLD(1 : ND, IND), 
     $                              JNEW(1 : ND, IND), 
     $                              SLNEW(1 : ND, IND),
     $                              POPNUM(1 : ND, LOW), 
     $                              POPNUM(1 : ND, NUP),
     $                              ARR(1 : ND, NUP, LOW), ARR(1 : ND, LOW, NUP),
     $                              ACR(1 : ND, NUP, LOW), ACR(1 : ND, LOW, NUP),
     $                              RBR(1 : ND, NUP, LOW))

         ENDDO

         DO I = 2, NLAST(1) - 1

            CALL PRINT_HYD_NLTE_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, T,
     $                              Z, Z, Z, Z, Z,
     $                              POPNUM(1 : ND, I), POPNUM(1 : ND,  NLAST(1)),
     $                              ARR(1 : ND, NLAST(1), I), ARR(1 : ND, I, NLAST(1)),
     $                              ACR(1 : ND, NLAST(1), I), ACR(1 : ND, I, NLAST(1)),
     $                              RBR(1 : ND, NLAST(1), I))

         ENDDO

         CALL RM_FILE(NTP_FILE, '-vf')

         IF (LAMBDA_ITER .EQ. 1) THEN

             CALL RM_FILE(NTW_FILE, '-vf')

             CALL OPEN_TO_APPEND(195, NTW_FILE)

         ENDIF

         DO IND = 1, LASTIND

            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            IF (LAMBDA_ITER .EQ. 1) WRITE(195, '(2(A10,2x),I2,2x,I2,2(2x,I3),2(2x,E15.7))') LEVEL(LOW),
     $                                                                                      LEVEL(NUP),
     $                                                                                      NOM(LOW),
     $                                                                                      NCHARG(LOW),
     $                                                                                      INT(WEIGHT(LOW)),
     $                                                                                      INT(WEIGHT(NUP)),
     $                                                                                      EINST(NUP, LOW),
     $                                                                                      XLAM

            CALL PRINT_NLTETRAPOP(POPNUM(1 : ND, LOW) * ENTOT(1 : ND),
     $                            POPNUM(1 : ND, NUP) * ENTOT(1 : ND),
     #                            DEPART_ZWAAN(1 : ND, LOW),
     $                            DEPART_ZWAAN(1 : ND, NUP))

         ENDDO

         IF (LAMBDA_ITER .EQ. 1) CLOSE(195)

      ENDIF

      IF (PRINT_LTE_ARR) THEN

         DO K = 1, NFEDGE(1); XJC_EDG(1 : ND) = XJC_EDG(1 : ND) + XJC_LTE(1 : ND, K) * FWEIGHT(K); ENDDO

         CALL PRINT_LTE_TRA('H MINUS..1', 'H I......1', 0,
     $                      POP_LTE(1 : ND, NFIRST(1)),
     $                      POP_LTE(1 : ND, NFIRST(1) + 1),
     $                      XJC_EDG(1 : ND),
     $                      ARR_LTE(1 : ND, NFIRST(1) + 1, 1),
     $                      ARR_LTE(1 : ND, 1, NFIRST(1) + 1))

         DO IND = 1, N_HI_LIN

            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            CALL PRINT_LTE_TRA(LEVEL(LOW), LEVEL(NUP), IND,
     $                         POP_LTE(1 : ND, LOW),
     $                         POP_LTE(1 : ND, NUP),
     $                         XJL_LTE(1 : ND, IND),
     $                         ARR_LTE(1 : ND, NUP, LOW),
     $                         ARR_LTE(1 : ND, LOW, NUP))

         ENDDO

         DO I = NFIRST(1) + 1, NLAST(1) - 1

            DO K = 1, NFEDGE(I); XJC_EDG(1 : ND) = XJC_EDG(1 : ND) + XJC_LTE(1 : ND, K) * FWEIGHT(K); ENDDO

            CALL PRINT_LTE_TRA(LEVEL(I), 'H II......', NLINE + I - 1,
     $                         POP_LTE(1 : ND, I),
     $                         POP_LTE(1 : ND, NLAST(1)),
     $                         XJC_EDG(1 : ND),
     $                         ARR_LTE(1 : ND, NLAST(1), I),
     $                         ARR_LTE(1 : ND, I, NLAST(1)))

         ENDDO

      ENDIF

      IF (ALLOCATED(ABXYZ_new)) DEALLOCATE(ABXYZ_new)

      IF (ALLOCATED(ARR_LTE)) DEALLOCATE(ARR_LTE)
      IF (ALLOCATED(XJL_LTE)) DEALLOCATE(XJL_LTE)
      IF (ALLOCATED(XJC_LTE)) DEALLOCATE(XJC_LTE)
      IF (ALLOCATED(POP_LTE)) DEALLOCATE(POP_LTE)
      IF (ALLOCATED(XJC_EDG)) DEALLOCATE(XJC_EDG)

      DEALLOCATE(ARR); DEALLOCATE(ACR); DEALLOCATE(RBR)

      DEALLOCATE(NOFILE)

      CLOSE(648); CLOSE(295); CLOSE(257)

!      INQUIRE(FILE = NRRM_FILE, EXIST = NRRM_FILE_EXISTS)
!      INQUIRE(FILE = NCRM_FILE, EXIST = NCRM_FILE_EXISTS)
!      INQUIRE(FILE = NTRM_FILE, EXIST = NTRM_FILE_EXISTS)

!      IF (NRRM_FILE_EXISTS) CALL REPLACE_PATTERN(NRRM_FILE, '0.0000000E+00', '      -      ', 13, 13)
!      IF (NCRM_FILE_EXISTS) CALL REPLACE_PATTERN(NCRM_FILE, '0.0000000E+00', '      -      ', 13, 13)
!      IF (NTRM_FILE_EXISTS) CALL REPLACE_PATTERN(NTRM_FILE, '0.0000000E+00', '      -      ', 13, 13)

      IF (LTE_RUN) STOP 'LTE RUN IS DONE'

      IF (.NOT. LTE_RUN) RETURN

      end subroutine

      subroutine print_nlte_lev(Level, LevPopLTE, LevPop, DepartZwaan)

      use file_operations
      use string_operations
      use common_block

      implicit none

      CHARACTER*10, INTENT(IN) ::           Level

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LevPop, LevPopLTE, DepartZwaan

      CHARACTER(:), ALLOCATABLE ::          FILE_NAME, LevName

      REAL*8 ::                             RAND

      integer ::                            file_unit

      integer ::                            di

      IF (Level(1:3) .NE. 'HEI') LevName = RM_CHAR(RM_CHAR(RM_CHAR(Level, ' '), '.'), '-')

      IF (Level(1:3) .EQ. 'HEI') LevName = RM_CHAR(Level, ' ')

      IF (LevName .EQ. 'HMINUS1')   LevName = 'HMINUS'
      IF (LevName .EQ. 'ELECTRONS') LevName = 'ELECTR'

      FILE_NAME = TRIM(ADJUSTL(NLTE_DIR_1//LevName))

      CALL RANDOM_NUMBER(RAND)

      FILE_UNIT = FLOOR(100D0 + RAND * 1000D0)

      if ((each_ali .and. lambda_iter .eq. 0) .or. .not. each_ali) then

          call rm_file(file_name, '-f')

          call open_to_append(file_unit, file_name)

          write(file_unit, '(3x,A,3x,A,8x,A,12x,A,14x,A,/)') 'li', 'di', 'lplte', 'lp', 'dz'

      else

          call open_to_append(file_unit, file_name)

      endif

      DO DI = 1, DPN

         WRITE(FILE_UNIT, '(I5,2x,I3,2x,E15.7,1x,E15.7,1x,E15.7)')
     $         LAMBDA_ITER, DI, LevPopLTE(DI), LevPop(DI), DepartZwaan(DI)

      ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_NLTE_LEV


      subroutine print_hyd_nlte_tra(ll, ul, idx, wvl, T, lo,
     $                              jold, sold, jnew, snew,
     $                              nl, nu, rul, rlu, cul, clu, rb)

      use file_operations
      use string_operations
      use common_block

      implicit none

      character*10, intent(in) ::           ul, ll

      integer, intent(in) ::                idx

      real*8, intent(in) ::                 wvl

      real*8, dimension(dpn), intent(in) :: lo, snew, jold, sold, jnew, rb

      real*8, dimension(dpn), intent(in) :: T, nl, nu

      real*8, dimension(dpn), intent(in) :: rul, rlu, cul, clu

      character(:), allocatable ::          file_name, tran_name, llev, ulev

      integer ::                            file_unit

      integer ::                            l

      character(len = 1000) ::              fmt_head, fmt_body

      fmt_head = '(A,9x,A,9x,A,4x,A,7x,A,8x,A,12x,A,18x,A,12x,A,8x,A,10(14x,A),/)'

      fmt_body = '(i3,2x,es15.7,2(2x,i5),2(3x,F9.2),2x,es15.7,2x,es23.15,4x,L,11(2x,es15.7))'

      llev = rm_char(rm_char(ll, ' '), '.')

      IF (llev .EQ. 'HMINUS1') llev = 'HMINUS'

      ulev = rm_char(rm_char(ul, ' '), '.')

      tran_name = llev//'_'//ulev

      file_name = trim(adjustl(nlte_dir_2//tran_name))

      file_unit = idx * 165

      if ((each_ali .and. lambda_iter .eq. 0) .or. .not. each_ali) then

           call rm_file(file_name, '-f')
            
           call open_to_append(file_unit, file_name)

           write(file_unit, fmt_head) 'wid', 'wvl', 'lit', 'hid', 'hei', 'tem', 'tau', 'llo',
     $                                'dam', 'jol', 'sol', 'jne', 'sne', 'nup', 'nlo',
     $                                'rul', 'rlu', 'cul', 'clu', 'rbr'

      else

           call open_to_append(file_unit, file_name)

      endif

      do l = 1, dpn

         write(file_unit, fmt_body) idx, wvl, lambda_iter, l, 
     $         height(l), T(l), tau_line(l, idx), lo(l), damp_line(l, idx),
     $         jold(l), sold(l), jnew(l), snew(l),
     $         nu(l), nl(l),
     $         rul(l), rlu(l),
     $         cul(l), clu(l),
     $         rb(l)

      enddo

      close(file_unit)

      end subroutine print_hyd_nlte_tra


      SUBROUTINE PRINT_NLTETRAPOP(PopLow, PopUp, DEPLOW, DEPUP)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      REAL*8, DIMENSION(DPN), INTENT(IN) :: PopLow, PopUp, DEPLOW, DEPUP

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      FILE_UNIT = 1653

      CALL OPEN_TO_APPEND(FILE_UNIT, NTP_FILE)

      DO DI = 1, DPN; WRITE(FILE_UNIT, '(E15.7,3(2x,E15.7))') PopLow(DI), PopUp(DI), DEPLOW(DI), DEPUP(DI); ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_NLTETRAPOP


      SUBROUTINE PRINT_LTE_LEV(LEV_NAM, LEV_NUM, LEV_POP)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           LEV_NAM

      INTEGER, INTENT(IN) ::                LEV_NUM

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LEV_POP

      CHARACTER(:), ALLOCATABLE ::          LEV_NAME, SEP_FILE_NAME

      INTEGER ::                            SEP_FILE_UNIT

      INTEGER ::                            DI
 
      LEV_NAME = RM_CHAR(RM_CHAR(LEV_NAM, ' '), '.')

      IF (LEV_NAME .EQ. 'HMINUS1')   LEV_NAME = 'HMINUS'
      IF (LEV_NAME .EQ. 'ELECTRONS') LEV_NAME = 'ELECTR'

      SEP_FILE_NAME = TRIM(ADJUSTL(LTE_DIR_1//LEV_NAME))

      SEP_FILE_UNIT = LEV_NUM * 13

      CALL OPEN_TO_APPEND(SEP_FILE_UNIT, SEP_FILE_NAME)

      WRITE(SEP_FILE_UNIT, '(1x,A,5x,A,2x,A,10x,A,/)') 'LEV', 'NUM', 'DI', 'POP'

      DO DI = 1, DPN; WRITE(SEP_FILE_UNIT, '(A6,2(2x,I3),3x,E15.7)') LEV_NAME, LEV_NUM, DI, LEV_POP(DI); ENDDO

      CLOSE(SEP_FILE_UNIT)

      END SUBROUTINE


      SUBROUTINE PRINT_LTE_TRA(LowLevel, UpLevel, LineInd,
     $                         PopLow, PopUp, TranInt,
     $                         RadUpLow, RadLowUp)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           UpLevel, LowLevel

      INTEGER, INTENT(IN) ::                LineInd

      REAL*8, DIMENSION(DPN), INTENT(IN) :: PopLow, PopUp

      REAL*8, DIMENSION(DPN), INTENT(IN) :: TranInt

      REAL*8, DIMENSION(DPN), INTENT(IN) :: RadUpLow, RadLowUp

      CHARACTER(:), ALLOCATABLE ::          FILE_NAME, TRAN_NAME, LowLev, UpLev

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      LowLev = RM_CHAR(RM_CHAR(LowLevel, ' '), '.')

      IF (LowLev .EQ. 'HMINUS1') LowLev = 'HMINUS'

      UpLev = RM_CHAR(RM_CHAR(UpLevel, ' '), '.')

      TRAN_NAME = LowLev//'_'//UpLev

      FILE_NAME = TRIM(ADJUSTL(LTE_DIR_2//TRAN_NAME))

      FILE_UNIT = LineInd * 165

      CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

      WRITE(FILE_UNIT, '(3x,A,6x,A,3x,A,2x,A,10x,A,2(14x,A),2(13x,A),/)') 'LL',
     $                                                                    'UL',
     $                                                                    'IND',
     $                                                                    'DI',
     $                                                                    'PU',
     $                                                                    'PL',
     $                                                                    'INT',
     $                                                                    'RUL',
     $                                                                    'RLU'

      DO DI = 1, DPN

         WRITE(FILE_UNIT, '(2(A6,2x),2(I3,1x),5(1x,E15.7))')
     $         LowLev, UpLev, LineInd, DI,
     $         PopUp(DI), PopLow(DI), TranInt(DI),
     $         RadUpLow(DI), RadLowUp(DI)

      ENDDO

      close(file_unit)

      end subroutine

      end module
