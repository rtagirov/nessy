      MODULE MOD_LINPOP

      CONTAINS

      SUBROUTINE LINPOP(T,RNE,ENTOT,ITNE,POPNUM,DEPART_ZWAAN,POP1,
     $   N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EN,ENDELTA,EINST,LEVEL,
     $   XLAMBDA,FWEIGHT,XJC,NF,XJL,WCHARM,XJCAPP,SCOLD,XJLAPP,
     $   DM,DB,XRED,XBLUE,V1,V2,V4,V5,VOLD,GAMMAL,EPSILON,
     $   TNEW,NOTEMP,NODM,IADR19,MAXADR,
     $   DELTAC,GAMMAR,IPRICC,IPRILC,MODHEAD,JOBNUM,IFRRA,ITORA,
     $   RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $   VELO,GRADI,VDOP,PHI,PWEIGHT,
     $   SCOLIND,SCNEIND,OPACIND,INDNUP,INDLOW,LASTIND,
     $   OPAC,SCNEW,DOPA,DETA,OPAL,DOPAL,DETAL,SIGMAKI,
     $   ND,LSRAT,CRATE,RRATE,RATCO,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,
     $   LINE,ALTESUM,ETAC,NFEDGE,EXPFAC,NOM,NATOM,KODAT,NFIRST,
     $   NLAST,WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,JOBMAX)

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
!     DB(N,N) = BROYDEN MATRIX AT CURRENT DEPTH POINT
!     DM(N,N) = ONLY NEEDED IF NO BROYDEN MATRIX FOUND AT START OF NEW MODEL
!     IF THERE IS A DB MATRIX, DM IS USED AS HELP-MATRIX FOR ALGORITHM.
!     VOLD(J) = VEKTOR OF K-1 ITERATION USED FOR K ITERATION
!     V1(J) = B-VEKTOR COMPUTED IN COMA()
!     V2-V5(J) = HELP-VEKTOR FOR BROYDEN ALGORITHMUS

      use MOD_ACOPY
      USE MOD_BFCROSS
      USE MOD_DBCLOSE
      USE MOD_DBLOAD
      USE MOD_DBOPEN
      USE MOD_DBSAVE
      use MOD_ERF_INF
      USE MOD_ERROR
      use MOD_FLGRID
      USE MOD_INV
      use MOD_ISRCHFGT
      use MOD_LTEPOP
      use MOD_PRIRAT
      use MOD_VADD
      use MOD_VADDM
      use MOD_VDIFF
      use MOD_VMD
      use MOD_VMF
      use MOD_VMT
      use MOD_VMV
      use MOD_VSUB
      USE MOD_XRUDI
      use MOD_CCORE
      use MOD_COMA
      use ABUNDANCES
      use MOD_LIOP

      USE FILE_OPERATIONS
      USE COMMON_BLOCK
      USE LINPOP_AUX

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION V1(NDIM+2),V2(NDIM+2),EN(NDIM+2),ENDELTA(NDIM+2),
     $          V4(NDIM+2),V5(NDIM+2),VOLD(NDIM+2),
     $          DB(NDIM+2,NDIM+2),DM(NDIM+2,NDIM+2)

      DIMENSION T(ND),RNE(ND),ENTOT(ND),ITNE(ND),TNEW(ND)
      DIMENSION XLAMBDA(NF),EXPFAC(NF)
      DIMENSION ELEVEL(NDIM),EION(NDIM),NCHARG(NDIM),NFEDGE(NDIM)
      DIMENSION EINST(NDIM,NDIM)

      DIMENSION CRATE(NDIM,NDIM),RRATE(NDIM,NDIM)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND),SCOLIND(LASTIND)
      DIMENSION XRED(LASTIND),XBLUE(LASTIND)
      DIMENSION SCOLD(NF,ND)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      REAL*8 SIGMAKI(NF), ALPHA(*), SEXPO(*)
      character*8 :: agaunt(NDIM)
      REAL*8 RATCO(*)
      real*8 rnel, rnellte

      real*8, allocatable :: ABXYZ_new(:)

      INTEGER NOM(N)

      LOGICAL NOTEMP, NODM, NEWRAP, STRONG_CONV, WEAK_CONV

      COMMON / COMNEGI / NEGINTL
      CHARACTER modhead*104
      CHARACTER LEVEL(NDIM)*10
      CHARACTER*4 KEYCOL(NDIM,NDIM)
      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)

      REAL*8 :: POPHIIL, POPHML, POPHIL

      real*8 PHI(*),PWEIGHT(*), WEIGHT(*)
      logical LINE(*), LBKG
      integer XLBKG1, XLBKG2
      character*10,dimension(*) :: MAINPRO,MAINLEV
      real*8,dimension(*) :: FWEIGHT, WCHARM, XJCAPP, XJLAPP
      integer, dimension(*) :: IWARN, KODAT
      real*8, dimension(*) :: ALTESUM, COCO, DETA, DETAL, DOPA, DOPAL
      real*8, dimension(*) :: ETA, ETAC, GRADI
      real*8, dimension(*) :: OPA, OPAC, OPACIND, OPAL, RADIUS, SCNEIND
      real*8, dimension(*) :: SCNEW, THOMSON, VELO

      REAL*8, DIMENSION(ND, LASTIND) :: XJL, JNEW
      REAL*8, DIMENSION(ND, NF) ::      XJC

      REAL*8, DIMENSION(:), ALLOCATABLE :: XJC_EDG

      REAL*8, DIMENSION(N) ::     ENLTE
      REAL*8, DIMENSION(ND, N) :: POPNUM, POPNUM_LTE, POP1
      REAL*8, DIMENSION(ND, N) :: DEPART_ZWAAN, DEPART_MENZEL

      REAL*8, DIMENSION(ND, LASTIND) :: SLOLD, SLNEW

      REAL*8, DIMENSION(1) :: ETAL1, OPAL1

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
      DEPART_MENZEL(1 : ND, 1 : N) = 0.0D0

      JNEW(1 : ND, 1 : LASTIND) = 0.0D0

      ElecConc(1 : ND) = 0.0D0
      ElecConcLTE(1 : ND) = 0.0D0

      ELEC_CONC_OLD(1 : ND) = 0.0D0

      ElecConcDep(1 : ND) = 0.0D0

      ALLOCATE(ARR(ND, N, N)); ARR(:, :, :) = 0.0D0
      ALLOCATE(ACR(ND, N, N)); ACR(:, :, :) = 0.0D0
      ALLOCATE(RBR(ND, N, N)); RBR(:, :, :) = 0.0D0

      Z(1 : ND) = 0.0D0
      ONE(1 : ND) =   1.0D0

C***  ENERGY EQUATION INCREASES THE RANK OF THE SYSTEM
      NPLUS1 = N + 1
      NPLUS2 = N + 2

      IF (NOTEMP) NRANK = NPLUS1

      IF (.NOT. NOTEMP) NRANK = NPLUS2

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
               CALL LIOP(EINST(NUP, LOW), WEIGHT(LOW), WEIGHT(NUP), LOW, NUP,
     $                   1, XLAM, ENTOT(L : L), POP1(L, 1 : N), RSTAR, OPAL1(1 : 1),
     $                   ETAL1(1 : 1), VDOP, N)

               IF (OPAL1(1) .LE. 0.0D0) THEN

                  SLOLD(L, IND) = 0.0D0

               ELSE

                  SLOLD(L, IND) = ETAL1(1) / OPAL1(1)

               ENDIF

            ENDIF

         ENDDO

      ENDDO

!     GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
!     SIGMAKI(K,LOW) IN CM**2

      CALL BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,NDIM,
     $             XLAMBDA,ALPHA,SEXPO,AGAUNT,NOM,WAVARR,SIGARR,NFDIM)

!     DETERMINE THE WEIGHT FUNCTION WCHARM FOR THE CONTINUUM, AND SCOLD AT ALL DEPTH POINTS
      CALL CCORE(WCHARM,NF,DELTAC,IPRICC,MODHEAD,JOBNUM,
     $           SCOLD,RADIUS,XLAMBDA,ND,T,RNE,POP1,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $           NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,SIGMAKI,
     $           WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,NFDIM)

!     GENERATE LINE FREQUENCY GRID AND INTEGRATION WEIGHTS
      CALL FLGRID(NFLDIM, NFL, PHI, PWEIGHT, DELTAX)
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

         CALL PRINT_HYD_NLTE_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, VDOP,
     $                           Z, Z, Z, Z, Z,
     $                           POP1(1 : ND, 1), POP1(1 : ND, 2),
     $                           Z, Z, Z, Z, Z)

         DO IND = 1, N_HI_LIN
 
            LOW = INDLOW(IND)
            NUP = INDNUP(IND)

            XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

            CALL PRINT_HYD_NLTE_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, VDOP,
     $                              LOO(1 : ND, IND), 
     $                              XJL(1 : ND, IND),
     $                              SLOLD(1 : ND, IND),
     $                              Z,
     $                              SLOLD(1 : ND, IND),
     $                              POP1(1 : ND, LOW), POP1(1 : ND, NUP),
     $                              Z, Z, Z, Z, Z)

         ENDDO

         DO I = 2, NLAST(1) - 1

            CALL PRINT_HYD_NLTE_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, VDOP,
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

      IF (.NOT. NOTEMP) EN(NPLUS2) = TL

!     IN CASE OF NON - ZERO LINE CORES CALCULATE BACKGROUND CONTINUUM SOURCE FUNCTION BY INTERPOLATION
      DO 11 IND = 1, LASTIND

         LOW = INDLOW(IND)
         NUP = INDNUP(IND)

         IF (EINST(LOW, NUP) .EQ. -2.D0) CYCLE

         WAVENUM = ELEVEL(NUP) - ELEVEL(LOW)

         CALL XRUDI(SCOLIND(IND),WAVENUM,SCOLD(1, L),XLAMBDA,1,NF,1)

   11 ENDDO

!     LOAD BROYDEN-MATRIX OF CURRENT DEPTH POINT
      IF (.NOT. NOFILE(L)) CALL DBLOAD(DB, L, NRANK)

      !***  BEGINNING OF THE BROYDEN ITERATION LOOP
      ITNE(L)=0
   10 ITNE(L)=ITNE(L)+1
      !***  UPDATE OF THE TEMPERATURE TL
      IF (ITNE(L) .GT. 1 .AND. .NOT. NOTEMP) TL=EN(NPLUS2)

      !***  PRE-CALCULATE EXPONENTIAL FACTORS FOR THE TEMPERATURE OF THE
      !***  CURRENT DEPTH POINT
      !***  THIS MUST BE REPEATED, WHEN THE TEMPERATURE HAS BEEN UPDATED
      IF (ITNE(L) .EQ. 1 .OR. .NOT. NOTEMP) THEN
      DO 25 K=1,NF
        WAVENUM=1.E8/XLAMBDA(K)
        EXPFAC(K)=EXP(-C1*WAVENUM/TL)
   25 CONTINUE
      ENDIF

      !***  CALCULATE LTE POP. NUMBERS
      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ_new,NFIRST,NLAST,NATOM)

!     CALCULATE NODM, NODM = TRUE => DM IS NOT NEEDED
      IF ((ITNE(L) .EQ. 1 .AND. NOFILE(L)) .OR. NEWRAP) THEN; NODM = .FALSE.; ELSE; NODM = .TRUE.; ENDIF

!     SETUP COEFFICIENT MATRICES
      POPHIIL = POPNUM(L, NLAST(1)) * ENTOT(L)
      POPHML  = POPNUM(L, 1)  *       ENTOT(L)
      POPHIL =  POPNUM(L, 2)  *       ENTOT(L)

!     LOO is declared in common_block.for along with its description

!      do iii = 1, LASTIND

!         print*, 'before:', iii, LOO(L, iii)

!      enddo

!      stop

      CALL COMA(CRATE,RRATE,RATCO,DM,N,NRANK,NDIM,V1,ABXYZ_new,
     $          ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,EION,WEIGHT,ALTESUM,
     $          XLAMBDA,FWEIGHT,XJC,NF, L, XJL(L, 1 : LASTIND), ND, XJLAPP,
     $          SLOLD(L, 1 : LASTIND), LASTIND, INDLOW,
     $          INDNUP,NOM,NATOM,KODAT,NFIRST,NLAST,PHI,PWEIGHT,DELTAX,XMAX,
     $          NFL,OPAC,SCNEW,DOPA,DETA,OPAL,SLNEW(L, 1 : LASTIND),
     $          DOPAL, DETAL, SIGMAKI,
     $          ETAC,NFEDGE,EXPFAC,SCOLIND,SCNEIND,OPACIND,NOTEMP,NODM,
     $          XRED,XBLUE,WCHARM,EN,RSTAR,SCOLD,XJCAPP,VDOP,COCO,KEYCOL,
     $          POPHIIL,POPHML, POPHIL, LOO(L, 1 : LASTIND), ITNE(L),
     $          LEVEL, JOBNUM, IRESTA)

!      do i = 1, 114

!         print*, 'after:', i, EN(i)

!      enddo

!      stop

      JNEW(L, 1 : LASTIND) = XJLAPP(1 : LASTIND)

      IF (NEWRAP) THEN

!       ALGEBRA OF ONE NEWTON ITERATION STEP
        CALL VMF(V2, EN, RATCO, NRANK, NRANK)

        CALL VSUB(V1, V2, NRANK) ! V1 = V1 - V2

        CALL INV(NRANK, NRANK, DM)

        CALL VMF(ENDELTA, V1, DM, NRANK, NRANK) ! ENDELTA = DM * V1

      ELSE

        !***  ALGEBRA OF ONE BROYDEN ITERATION STEP
        !***  calculate the resulting vector V4 using the current populations EN

        CALL VMF(V4, EN, RATCO, NRANK, NRANK)

        !*** calculate the error vector V1 = V1 - V4

        CALL VSUB (V1,V4,NRANK)

        IF (ITNE(L) .EQ. 1) THEN

           IF (NOFILE(L)) THEN

!           IF FORT.19 NOT FOUND AND THE FIRST ITERATION STEP, DB IS DM^T
            CALL INV (NRANK,NRANK,DM)
            CALL ACOPY (DB,DM,NRANK,NRANK)

           ENDIF

        ELSE

C***  ALGEBRA FOR THE BROYDEN FORMULA
c       delta_y = VOLD = VOLD-V1
            CALL VSUB(VOLD,V1,NRANK)

c       V2 = delta_y * B
            CALL VMF(V2,VOLD,DB,NRANK,NRANK)

c       V3 = V2*ENDELTA; ENDELTA is the difference between pops = delta_x
c       V3 the denominater (a number)
            CALL VMV(V3,V2,ENDELTA,NRANK,NRANK)

            V5(:NRANK)=ENDELTA(:NRANK)

c       ENDELTA = B*delta_x = DB*V5
            CALL VMT(ENDELTA,DB,V5,NRANK,NRANK)

c       V5 = delta_x - delta_y*B = V5-V2
            CALL VSUB (V5,V2,NRANK)   ! V5=V5-V2

c       DM = B*delta_x * V5
c       numerator (matrix) DM = dyadic product ENDELTA * V5
            CALL VMD (DM,ENDELTA,V5,NRANK,NRANK)

c       matrix DM (numerator) divided by V3
            CALL VDIFF (DM,V3,NRANK,NRANK)

c       new Broyden matrix DB = DB+DM
            CALL VADDM (DB,DM,NRANK,NRANK)

         ENDIF

         VOLD(:NRANK) = V1(:NRANK)

C      calculate the new population-correction vector ENDELTA = V1*DB
         CALL VMF(ENDELTA, V1, DB, NRANK, NRANK)

      ENDIF

!     ALGEBRA COMMON TO BOTH APPROACHES
!     UPDATE THE POPULATIONS EN = EN + ENDELTA
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
      IF (NEWRAP) CALL ACOPY(DB, DM, NRANK, NRANK)

      CALL DBSAVE(DB, L, NRANK)

!     UPDATING ARRAY POPNUM AND CALCULATING THE DEPARTURE COEFFICIENTS:
!     RINAT TAGIROV:
!     THE MOST ACCURATE WAY TO DO THIS IS TO UPDATE THE LTE POPULATIONS HERE BECAUSE UPON THE LAST
!     CALL OF THE LTEPOP PROCEDURE THE ELECTRON CONCENTRATION HAD NOT
!     YET CONVERGED, I.E. IF WE DON'T DO THAT HERE THE LTE POPULATIONS WILL
!     CORRESPOND TO THE PENULTIMATE BROYDEN/NEWTON ITERATION
!*******************************************************************************************

      IF (CONST_ELEC) THEN ! PRE-SET ELECTRON CONCENTRATION

          ELEC_CONC =       READ_ATM_MOD(ATM_MOD_FILE, '3')
          HEAVY_ELEM_CONC = READ_ATM_MOD(ATM_MOD_FILE, '4')

          EN(NPLUS1) = ELEC_CONC(L) / HEAVY_ELEM_CONC(L)

      ENDIF

      ENE = EN(NPLUS1) * ENTOT(L)

      CALL LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ_new,NFIRST,NLAST,NATOM)
!*******************************************************************************************

      POPNUM(L, 1 : N) =       EN(1 : N)
      POPNUM_LTE(L, 1 : N) =   ENLTE(1 : N)
      DEPART_ZWAAN(L, 1 : N) = EN(1 : N) / ENLTE(1 : N)

!      DEPART_MENZEL(L, 2 : NLAST(1) - 1) = DEPART_ZWAAN(L, 2 : NLAST(1) - 1) / DEPART_ZWAAN(L, NLAST(1))

C***  UPDATING THE ELECTRON DENSITY
      RNE(L) = EN(NPLUS1)

C***  UPDATING THE ELECTRON TEMPERATURE

      IF (.NOT. NOTEMP) THEN; TNEW(L) = EN(NPLUS2); ELSE; TNEW(L) = T(L); ENDIF

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
     $CALL       PRIRAT (ITNE(L),N,LEVEL,NDIM,L,CRATE,RRATE,RATCO,EN,
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

         CALL PRINT_HYD_NLTE_TRA('H MINUS..1', 'H I......1', 0, 0.0D0, VDOP,
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

            CALL PRINT_HYD_NLTE_TRA(LEVEL(LOW), LEVEL(NUP), IND, XLAM, VDOP,
     $                              LOO(1 : ND, IND), 
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

            CALL PRINT_HYD_NLTE_TRA(LEVEL(I), 'H II......', NLINE + I - 1, 0.0D0, VDOP,
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

!        RINAT TAGIROV: WRITING THE NEW ELECTRON CONCENTRATION TO THE UPDATED ATMOSPHERE MODEL FILE FAL_VD.UPD
         CALL UPDATE_ATM_MOD(RNE(1 : ND) * ENTOT(1 : ND))

      ENDIF

!      WRITE(*, '(/,A)') '******** REFINEMENT PROTONS **********'

!      DO L = 1, ND

!         IF (L .LE. 49 .AND. L .GE. 42) THEN

!            ONE_PRO = 0.0D0
!            ONE_PRO = POPNUM(L, 2) * ACR(L, 2, 12) - POPNUM(L, 12) * ACR(L, 12, 2)
!            ONE_PRO = ONE_PRO - POPNUM(L, 12) * RBR(L, 12, 2)
            
!            TWO_PRO = -1.0D0 * POPNUM(L, 12) * RBR(L, 12, 3)
!            TWO_PRO = TWO_PRO + POPNUM(L, 3) * ACR(L, 3, 12) - POPNUM(L, 12) * ACR(L, 12, 3)

!            PRO_LEV = 0.0D0

!            DO I = 4, 8; PRO_LEV = PRO_LEV + POPNUM(L, 12) * RBR(L, 12, I); ENDDO

!            DO I = 4, 11
!
!               PRO_LEV = PRO_LEV + POPNUM(L, 12) * RBR(L, 12, I) +
!     $                             POPNUM(L, 12) * ACR(L, 12, I) -
!     $                             POPNUM(L, I) *  ACR(L, I, 12)

!            ENDDO

!            WRITE(*, '(I3,2(1x,E15.7))') L, ONE_PRO, PRO_LEV * (1.0D0 - TWO_PRO / PRO_LEV)
!
!         ENDIF
!
!      ENDDO

!      WRITE(*, '(/,A)') '******** REFINEMENT 2ND **********'
!
!      DO L = 1, ND
!
!         IF (L .LE. 49 .AND. L .GE. 42) THEN
!
!            LEV_TWO = POPNUM(L, 4) * RBR(L, 4, 3) + POPNUM(L, 5) * RBR(L, 5, 3)
!            
!            TWO_PRO = -1.0D0 * POPNUM(L, 3) * RBR(L, 12, 3)
!
!            TWO_ONE = POPNUM(L, 3) * ACR(L, 3, 2) - POPNUM(L, 2) * ACR(L, 2, 3)
!
!            WRITE(*, '(I3,5(1x,E15.7))') L, LEV_TWO - TWO_PRO - TWO_ONE, TWO_ONE, LEV_TWO - TWO_PRO, LEV_TWO, TWO_PRO
!
!         ENDIF
!
!      ENDDO

!      WRITE(*, '(A,/)') '**********************************'

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

      END SUBROUTINE LINPOP

      END MODULE MOD_LINPOP
