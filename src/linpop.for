      module mod_linpop

      contains

      subroutine LINPOP(T,
     $                  RNE,
     $                  ENTOT,
     $                  ITNE,
     $                  POPNUM,
     $                  DEPART,
     $                  pop1_full,
     $                  N,
     $                  WEIGHT,
     $                  NCHARG,
     $                  EION,
     $                  ELEVEL,
     $                  EINST,
     $                  LEVEL,
     $                  XLAMBDA,
     $                  FWEIGHT,
     $                  XJC,
     $                  NF,
     $                  XJL,
     $                  WCHARM,
     $                  EPSILON,
     $                  MODHEAD,
     $                  JOBNUM,
     $                  IFRRA,
     $                  ITORA,
     $                  RADIUS,
     $                  RSTAR,
     $                  IWARN,
     $                  MAINPRO,
     $                  MAINLEV,
     $                  VDOP,
     $                  INDNUP,
     $                  INDLOW,
     $                  LASTIND,
     $                  ND,
     $                  LSRAT,
     $                  ALPHA,
     $                  SEXPO,
     $                  AGAUNT,
     $                  COCO,
     $                  KEYCOL,
     $                  ALTESUM,
     $                  NOM,
     $                  NATOM,
     $                  KODAT,
     $                  levatnum,
     $                  NFIRST,
     $                  NLAST,
     $                  WAVARR,
     $                  SIGARR,
     $                  LBKG,
     $                  XLBKG1,
     $                  XLBKG2,
     $                  JOBMAX,
     $                  N_full,
     $                  weight_full,
     $                  ncharg_full,
     $                  eion_full,
     $                  elevel_full,
     $                  einst_full,
     $                  level_full,
     $                  alpha_full,
     $                  sexpo_full,
     $                  agaunt_full,
     $                  coco_full,
     $                  keycol_full,
     $                  altesum_full,
     $                  nom_full,
     $                  natom_full,
     $                  kodat_full,
     $                  nfirst_full,
     $                  nlast_full)

!     this procedure follows Koesterke et al. 1992 A&A 255, 490

!     CALCULATION OF NEW NLTE POPULATION NUMBERS (ARRAY POPNUM)
!     RADIATIVE RATES ARE CALCULATED WITH THE SCHARMER RADIATION FIELD
!     THE (HENCE NON-LINEAR) RATE EQUATIONS ARE SOLVED BY LINEARIZATION
!     pop1_full = OLD LTE AND NLTE POPULATION NUMBERS
!     pop1 =      OLD NLTE POPULATION NUMBERS
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
!     V2 - V5(J) = HELP-VEKTOR FOR BROYDEN ALGORITHM

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
      use MOD_LIOP

      use file_operations
      use common_block
      use matoper
      use broyden
      use vardatom_lte

      implicit real*8(A - H, O - Z)

      integer,       intent(in)                                   :: lastind, N

      real*8,        intent(in),     dimension(N, N)              :: einst

      real*8,        intent(in),     dimension(N)                 :: elevel, eion, alpha, sexpo, weight

      integer,       intent(in),     dimension(N)                 :: ncharg, nom, levatnum

      integer,       intent(in),     dimension(natom)             :: nfirst, nlast, kodat

      integer,       intent(in),     dimension(lastind)           :: indnup, indlow

      character*8,   intent(in),     dimension(N)                 :: agaunt

      character*4,   intent(in),     dimension(N, N)              :: keycol

      character*10,  intent(in),     dimension(N)                 :: level

      real*8,        intent(in),     dimension(N, N, 4)           :: coco
      real*8,        intent(in),     dimension(4, N)              :: altesum

      real*8,        intent(in),     dimension(N_full, N_full)    :: einst_full

      real*8,        intent(in),     dimension(N_full)            :: elevel_full, eion_full

      real*8,        intent(in),     dimension(N_full)            :: alpha_full, sexpo_full, weight_full

      integer,       intent(in),     dimension(N_full)            :: ncharg_full, nom_full

      integer,       intent(in),     dimension(natom_full)        :: nfirst_full, nlast_full, kodat_full

      character*8,   intent(in),     dimension(N_full)            :: agaunt_full

      character*4,   intent(in),     dimension(N_full, N_full)    :: keycol_full

      character*10,  intent(in),     dimension(N_full)            :: level_full

      real*8,        intent(in),     dimension(N_full, N_full, 4) :: coco_full
      real*8,        intent(in),     dimension(4, N_full)         :: altesum_full

      real*8,        intent(in),     dimension(N_full, NF)        :: WAVARR, SIGARR

      real*8,        intent(in),     dimension(ND, N_full)        :: pop1_full

      real*8,        intent(in),     dimension(ND)                :: T, ENTOT, radius

      real*8,        intent(in),     dimension(NF)                :: xlambda, fweight
      real*8,        intent(in),     dimension(ND, NF)            :: WCHARM
      REAL*8,        intent(in),     dimension(ND, NF)            :: XJC

      character*10,  intent(in),     dimension(ND)                :: MAINPRO, MAINLEV
      integer,       intent(in),     dimension(ND)                :: IWARN
 
      character*104, intent(in)                                   :: modhead

      logical,       intent(in)                                   :: LBKG

      integer,       intent(in)                                   :: XLBKG1, XLBKG2
      integer,       intent(in)                                   :: JOBNUM

      integer,       intent(out),    dimension(ND)                :: itne

      real*8,        intent(out),    dimension(ND, N_full)        :: POPNUM

      real*8,        intent(inout),  dimension(ND, lastind)       :: XJL

      real*8,        intent(inout),  dimension(ND)                :: RNE

      integer,                       dimension(N_full)            :: nfedge

      real*8,                        dimension(ND, N)             :: POPNUM_nlte

      real*8,                        dimension(N + 1)             :: V1, V2, ENDELTA, V4, V5, VOLD
      real*8,                        dimension(N + 1, N + 1)      :: DB, DM

      real*8,                        dimension(N)                 :: ENLTE

      real*8,                        dimension(N_full)            :: ENLTE_full, en_full

      real*8,                        dimension(ND, N_full)        :: POPNUM_LTE

      real*8,                        dimension(N, N)              :: CRATE, RRATE

      real*8,                        dimension(N + 1)             :: en

      real*8,                        dimension(NF)                :: EXPFAC

      real*8,                        dimension(NF, ND)            :: SCOLD

      real*8,                        dimension(NF, N_full)        :: SIGMAKI

      real*8,                        dimension(N + 1, N + 1)      :: RATCO

      real*8,                        dimension(ND, lastind)       :: JNEW

      real*8,                        dimension(lastind)           :: XJLAPP

      real*8,                        dimension(ND, N_full)        :: DEPART

      real*8,                        dimension(ND, N)             :: pop1

      real*8,                        dimension(ND, LASTIND)       :: SLOLD, SLNEW

      real*8,                        dimension(ND)                :: ELEC_CONC_OLD, ElecConc, ElecConcLTE, ElecConcDep

      real*8,                        dimension(ND)                :: ONE, Z

      real*8,        allocatable,    dimension(:)                 :: XJC_EDG

      character(:),  allocatable                                  :: CONV_FILE

      integer                                                     :: IND, LOW, NUP, N_HI_LEV, N_HI_LIN

      real*8                                                      :: epsdn

      real*8                                                      :: POPHIIL, POPHML, POPHIL

      real*8                                                      :: opalind, etalind

      logical                                                     :: NODM, NEWRAP, STRONG_CONV, WEAK_CONV

      logical                                                     :: NRRM_FILE_EXISTS, NCRM_FILE_EXISTS

      logical                                                     :: NTRM_FILE_EXISTS, CONV_FILE_EXISTS

      logical                                                     :: PRINT_LTE_ARR

      integer :: iii

!     C1 = H * C / K (CM * KELVIN)
      DATA C1 /1.4388D0/

      CONV_FILE = CONV_DIR//'ALL'

      PRINT_LTE_ARR = .FALSE.

      DEPART(1 : ND, 1 : N_full) = 0.0D0

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

      SLOLD(1 : ND, 1 : LASTIND) = 0.0D0

!     REMOVE NEGATIVE LINE INTENSITIES
      NEGINTL = 0

      DO L = 1, ND

         do j = 1, N; pop1(L, j) = pop1_full(L, idx_orig(j)); enddo

         DO IND = 1, LASTIND

            LOW = INDLOW(IND)

            NUP = INDNUP(IND)

            IF (EINST(LOW, NUP) .NE. -2.0D0) THEN

               IF (XJL(L, IND) .LT. 0.0D0) THEN

                  XJL(L, IND) = ABS(XJL(L, IND))
 
                  NEGINTL = NEGINTL + 1

               ENDIF

               XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

!              calculate the opacity (opalind) and emissivity (etalind) for a given line and depth point
               CALL LIOP_SBE(EINST(NUP, LOW), WEIGHT(LOW), WEIGHT(NUP), LOW, NUP,
     $                       XLAM, ENTOT(L), pop1(L, 1 : N), RSTAR, opalind, etalind, VDOP, N)

               IF (opalind .LE. 0.0D0) THEN

                  SLOLD(L, IND) = 0.0D0

               ELSE

                  SLOLD(L, IND) = etalind / opalind

               ENDIF

            ENDIF

         ENDDO

      ENDDO

!     GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES SIGMAKI(K, LOW) IN CM**2
      CALL BFCROSS(SIGMAKI,NF,N_full,NCHARG_full,ELEVEL_full,EION_full,EINST_full,
     $             XLAMBDA(1 : NF),ALPHA_full,SEXPO_full,AGAUNT_full,NOM_full,
     $             WAVARR(:, 1 : NF),SIGARR(:, 1 : NF))

!     DETERMINE SCOLD AT ALL DEPTH POINTS
      CALL CCORE(NF,MODHEAD,JOBNUM,SCOLD,RADIUS,XLAMBDA,ND,T,RNE,pop1_full,ENTOT,RSTAR,
     $           IWARN,MAINPRO,MAINLEV,NOM_full,
     $           N_full,LEVEL_full,NCHARG_full,WEIGHT_full,ELEVEL_full,EION_full,EINST_full,SIGMAKI,
     $           WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2)

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

         DO J = 1, N_full; CALL PRINT_LEV(LEVEL_full(J), pop1_full(1 : ND, J), pop1_full(1 : ND, J), ONE); ENDDO

         CALL PRINT_LEV('ELECTRONS ', RNE(1 : ND), RNE(1 : ND), ONE)

         CALL PRINT_HYD_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, T,
     $                      Z, Z, Z, Z, Z,
     $                      pop1(1 : ND, 1),
     $                      pop1(1 : ND, 2),
     $                      Z, Z, Z, Z, Z)

         DO IND = 1, N_HI_LIN
 
            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            CALL PRINT_HYD_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, T,
     $                         LLO(1 : ND, IND),
     $                         XJL(1 : ND, IND),
     $                         SLOLD(1 : ND, IND),Z,
     $                         SLOLD(1 : ND, IND),
     $                         pop1(1 : ND, LOW),
     $                         pop1(1 : ND, NUP),
     $                         Z, Z, Z, Z, Z)

         ENDDO

         DO I = 2, NLAST(1) - 1

            CALL PRINT_HYD_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, T,
     $                         Z, Z, Z, Z, Z,
     $                         pop1(1 : ND, I),
     $                         pop1(1 : ND, NLAST(1)),
     $                         Z, Z, Z, Z, Z)

         ENDDO

         LAMBDA_ITER = 1

         IF (JOBMAX .LE. 5.) STOP 'LTE IS DONE.'

      ELSEIF (LTE_RUN) THEN

         CALL MKDIR(LTE_DIR_1); CALL CLEAN_DIR(LTE_DIR_1)
         CALL MKDIR(LTE_DIR_2); CALL CLEAN_DIR(LTE_DIR_2)

         ALLOCATE(XJL_LTE(ND, LASTIND)); XJL_LTE(1 : ND, 1 : LASTIND) = XJL(1 : ND,  1 : LASTIND)
         ALLOCATE(XJC_LTE(ND, NF));      XJC_LTE(1 : ND, 1 : NF) =      XJC(1 : ND,  1 : NF)
         ALLOCATE(POP_LTE(ND, N));       POP_LTE(1 : ND, 1 : N) =       pop1(1 : ND, 1 : N)

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

      ELEC_CONC_OLD(L) = SUM(NCHARG_full(1 : N_full) * pop1_full(L, 1 : N_full))

!     CALCULATE FREQUENCY INDICES OF IONIZATION EDGES
      DO NA = 1, NATOM_full

         DO LOW = NFIRST_full(NA), NLAST_full(NA) - 1

            EDGE = EION_full(LOW) - ELEVEL_full(LOW)

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

      EN(1 : N) = pop1(L, 1 : N)

      EN(NPLUS1) = RNE(L)

      en_full(1 : N_full) = pop1_full(L, 1 : N_full)

      TL = T(L)

!     LOAD BROYDEN-MATRIX OF CURRENT DEPTH POINT
      IF (.NOT. NOFILE(L)) CALL DBLOAD(DB, L, NRANK)

      !BEGINNING OF THE BROYDEN ITERATION LOOP
      ITNE(L) = 0
   10 ITNE(L) = ITNE(L) + 1

!     PRE-CALCULATE EXPONENTIAL FACTORS FOR THE TEMPERATURE OF THE CURRENT DEPTH POINT
!     THIS MUST BE REPEATED, WHEN THE TEMPERATURE HAS BEEN UPDATED
      if (ITNE(L) .eq. 1) then

          do K = 1, NF

             WAVENUM = 1.E8 / XLAMBDA(K)

             EXPFAC(K) = EXP(-C1 * WAVENUM / TL)

          enddo

      endif

!     CALCULATE LTE POPULATION NUMBERS
      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZn_nlte(1 : NATOM, L),NFIRST,NLAST,NATOM)

!     CALCULATE NODM, NODM = TRUE => DM IS NOT NEEDED
      IF ((ITNE(L) .EQ. 1 .AND. NOFILE(L)) .OR. NEWRAP) THEN; NODM = .FALSE.; ELSE; NODM = .TRUE.; ENDIF

!     SETUP COEFFICIENT MATRICES
      POPHIIL = POPNUM(L, NLAST_full(1)) * ENTOT(L)
      POPHML  = POPNUM(L, 1) *             ENTOT(L)
      POPHIL =  POPNUM(L, 2) *             ENTOT(L)

      CALL COMA(CRATE,
     $          RRATE,
     $          RATCO,
     $          DM,
     $          N,
     $          NRANK,
     $          V1,
     $          ABXYZn_nlte(1 : NATOM, L),
     $          ENLTE,
     $          TL,
     $          ENE,
     $          NCHARG,
     $          ELEVEL,
     $          EINST,
     $          EION,
     $          WEIGHT,
     $          ALTESUM,
     $          XLAMBDA,
     $          FWEIGHT,
     $          XJC,
     $          NF,
     $          L,
     $          XJL(L, 1 : LASTIND),
     $          ND,
     $          XJLAPP,
     $          SLOLD(L, 1 : LASTIND),
     $          LASTIND,
     $          INDLOW,
     $          INDNUP,
     $          NOM,
     $          NATOM,
     $          KODAT,
     $          levatnum,
     $          NFIRST,
     $          NLAST,
     $          SLNEW(L, 1 : LASTIND),
     $          SIGMAKI,
     $          NFEDGE,
     $          EXPFAC,
     $          NODM,
     $          WCHARM,
     $          EN,
     $          RSTAR,
     $          SCOLD,
     $          VDOP,
     $          COCO,
     $          KEYCOL,
     $          POPHIIL,
     $          POPHML,
     $          POPHIL,
     $          LLO(L, 1 : LASTIND), ! LLO is declared in comblock.for
     $          ITNE(L),
     $          LEVEL,
     $          JOBNUM,
     $          IRESTA,
     $          N_full, ncharg_full, weight_full, elevel_full, eion_full, en_full, nom_full)

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

!     updating array popnum and calculating the departure coefficients:
!     Rinat Tagirov:
!     the most accurate way to do this is to update the lte populations here
!     because upon the last LTEPOP call the electron concentration has not
!     yet converged, i.e. if we don't do that here the LTE populations will
!     correspond to the penultimate broyden/newton iteration
!*******************************************************************************************

      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N_full,
     $            ENLTE_full,
     $            TL,ENE,
     $            WEIGHT_full,
     $            NCHARG_full,
     $            EION_full,
     $            ELEVEL_full,
     $            NOM_full,
     $            ABXYZn(1 : NATOM_full, L),
     $            NFIRST_full,
     $            NLAST_full,
     $            NATOM_full)

      do i = 1, N; en_full(idx_orig(i)) = EN(i); enddo

      if (natom_lte /= 0) then

          do i = 1, N_full; if (.not. nlte_lev(i)) en_full(i) = enlte_full(i); enddo

      endif

!*******************************************************************************************

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

      POPNUM_NLTE(L, 1 : N) = EN(1 : N)

      POPNUM_LTE(L, 1 : N_full) = ENLTE_full(1 : N_full)
      POPNUM(L,     1 : N_full) = ENLTE_full(1 : N_full)

!     the populations of levels treated in NLTE are replaced with their NLTE values
      do i = 1, N; POPNUM(L, idx_orig(i)) = EN(i); enddo

      DEPART(L, 1 : N_full) = POPNUM(L, 1 : N_full) / POPNUM_LTE(L, 1 : N_full)

!     UPDATING THE ELECTRON DENSITY
      RNE(L) = EN(NPLUS1)

      ElecConc(L) =    SUM(NCHARG_full(1 : N_full) * POPNUM(L,     1 : N_full))
      ElecConcLTE(L) = SUM(NCHARG_full(1 : N_full) * POPNUM_LTE(L, 1 : N_full))

      ElecConcDep(L) = ElecConc(L) / ElecConcLTE(L)

!     PRINTOUT OF RATE COEFFICIENTS
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

  100 CONTINUE

!     CLOSE FORT19 FOR BROYDEN-MATRIX
      CALL DBCLOSE(ND, NRANK, DB)

!     CALCULATING CORMAX FOR ELECTRONS
      CORMAX_ELEC = MAXVAL(DABS(1.0D0 - ELEC_CONC_OLD(3 : ND - 1) / ElecConc(3 : ND - 1)))

      IF (.NOT. LTE_RUN) THEN

         DO J = 1, N_full; CALL PRINT_LEV(LEVEL_full(J), POPNUM_LTE(1 : ND, J),
     $                                    POPNUM(1 : ND, J), DEPART(1 : ND, J)); ENDDO

         CALL PRINT_LEV('ELECTRONS ', ElecConcLTE(1 : ND), ElecConc(1 : ND), ElecConcDep(1 : ND))

         CALL PRINT_HYD_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, T,
     $                      Z, Z, Z, Z, Z,
     $                      POPNUM_nlte(1 : ND, NFIRST(1)),
     $                      POPNUM_nlte(1 : ND, NFIRST(1) + 1),
     $                      ARR(1 : ND, NFIRST(1) + 1, NFIRST(1)),
     $                      ARR(1 : ND, NFIRST(1),     NFIRST(1) + 1),
     $                      ACR(1 : ND, NFIRST(1) + 1, NFIRST(1)),
     $                      ACR(1 : ND, NFIRST(1),     NFIRST(1) + 1),
     $                      RBR(1 : ND, NFIRST(1) + 1, NFIRST(1)))

         DO IND = 1, N_HI_LIN
 
            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            CALL PRINT_HYD_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, T,
     $                         LLO(1 : ND, IND),
     $                         XJL(1 : ND, IND),
     $                         SLOLD(1 : ND, IND), 
     $                         JNEW(1 : ND, IND),
     $                         SLNEW(1 : ND, IND),
     $                         POPNUM_nlte(1 : ND, LOW),
     $                         POPNUM_nlte(1 : ND, NUP),
     $                         ARR(1 : ND, NUP, LOW),
     $                         ARR(1 : ND, LOW, NUP),
     $                         ACR(1 : ND, NUP, LOW),
     $                         ACR(1 : ND, LOW, NUP),
     $                         RBR(1 : ND, NUP, LOW))

         ENDDO

         DO I = 2, NLAST(1) - 1

            CALL PRINT_HYD_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, T,
     $                         Z, Z, Z, Z, Z,
     $                         POPNUM_nlte(1 : ND, I),
     $                         POPNUM_nlte(1 : ND, NLAST(1)),
     $                         ARR(1 : ND, NLAST(1), I),
     $                         ARR(1 : ND, I, NLAST(1)),
     $                         ACR(1 : ND, NLAST(1), I),
     $                         ACR(1 : ND, I, NLAST(1)),
     $                         RBR(1 : ND, NLAST(1), I))

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

            CALL PRINT_NLTETRAPOP(POPNUM_nlte(1 : ND, LOW) * ENTOT(1 : ND),
     $                            POPNUM_nlte(1 : ND, NUP) * ENTOT(1 : ND),
     #                            DEPART(1 : ND, LOW),
     $                            DEPART(1 : ND, NUP))

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

      subroutine print_lev(Level, LevPopLTE, LevPop, Depart)

      use file_operations
      use string_operations
      use common_block

      implicit none

      CHARACTER*10, INTENT(IN) ::           Level

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LevPop, LevPopLTE, Depart

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
     $         LAMBDA_ITER, DI, LevPopLTE(DI), LevPop(DI), Depart(DI)

      ENDDO

      CLOSE(FILE_UNIT)

      end subroutine


      subroutine print_hyd_tra(ll, ul, idx, wvl, T, lo,
     $                         jold, sold, jnew, snew,
     $                         nl, nu, rul, rlu, cul, clu, rb)

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

      end subroutine


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
