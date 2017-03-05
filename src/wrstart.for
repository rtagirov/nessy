      MODULE MOD_WRSTART

      CONTAINS

      SUBROUTINE WRSTART

      use MOD_DATOM_M
      use MOD_DECON
      use MOD_DECSTAR_M
      use MOD_FGRID
      use MOD_GEOMESH
      use MOD_GRADIFF
      use MOD_GREYM
      use MOD_PRIBLA
      use MOD_PRIDAT
      use MOD_PRIGH
      use MOD_PRIMOD
      use MOD_WRITPOP
      use MOD_WRVEL
      use MOD_JSTART
      use MOD_PRICOMP
      use MOD_REBLANK
      use MOD_WRITMOD
      use MOD_TICTOC
      use MOD_ERROR
      use MOD_chemeq      
      use ABUNDANCES
      USE COMMON_BLOCK
      USE FILE_OPERATIONS
      USE VARDATOM
      USE VARHMINUS

!     THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
!     CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
!     IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
!     AND THE FREQUENCY GRID (FILE FGRID)
!     PRESENT VERSION: MODEL ATMOSPHERE OF HELIUM (CODE NR. "1") WITH
!                                          HYDROGEN         "2"
!     FOR IMPLEMENTATION OF ADDITIONAL ELEMENTS:
!     MODIFY SUBROUTINES  "DATOM", "DECSTAR"
!     INSERT CORRESPONDING ATOMIC DATA INTO SUBR. "COLLI", "PHOTOCS"

      IMPLICIT REAL*8(A - H, O - Z)

      INTEGER :: IPDIM, NBDIM
      parameter(IPDIM = 25, NBDIM = 99)

      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM, NBDIM), ABSEVT(IPDIM, NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM, NFDIM), ABSFAC(NDDIM, NFDIM)
      COMMON /VELPAR/  VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE
      COMMON /COMTEFF/ TEFF, TMIN, TMODIFY, SPHERIC
      COMMON /COMLBKG/ LBKG, XLBKG1, XLBKG2

!     LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
!     XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF

      INTEGER XLBKG1, XLBKG2
      LOGICAL LBKG
      LOGICAL TTABLE, TPLOT, SPHERIC, FAL

      CHARACTER MODHEAD*104

      CHARACTER   NAME*10, fstring*24
      integer timer

      real*8 ATMEAN, AMU

      integer NA
      REAL*8 CSARR(5000,4)

      REAL*8, DIMENSION(:), ALLOCATABLE :: ELEC_CONC, HEAVY_ELEM_CONC

      REAL*8, DIMENSION(:), ALLOCATABLE :: VELO_NE, VELO_E

      CHARACTER(:), ALLOCATABLE :: AMF

      REAL*8 :: H

      DATA AMU /1.660531d-24/

      call FDATE(fstring)
      call TIC(timer)

!     READ ATOMIC DATA FROM FILE DATOM
      CALL DATOM_M(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $             EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $             INDNUP,INDLOW,LASTIND,NATOM,
     $             ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $             NLAST,WAVARR,SIGARR,NFDIM)

!     DECODING INPUT DATA
      CALL DECSTAR_M(MODHEAD,FM,RSTAR,VDOP,TTABLE,FAL,LBKG,XLBKG1,XLBKG2,
     $               TPLOT,NATOM,ABXYZ,KODAT,IDAT,LBLANK,ATMEAN, AMU)

!     if PRINT DATOM option in CARDS is set, printout the atomic data
      IF (IDAT.EQ.1)
     $CALL PRIDAT(N,LEVEL,NCHARG, WEIGHT,ELEVEL,EION,EINST,
     $            KODAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $            NATOM,ELEMENT,NOM,ABXYZ,ATMASS)

      !***  PRINTOUT OF THE CHEMICAL COMPOSITION
      CALL PRICOMP(N,EINST,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     $             STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     $             INDLOW,INDNUP)

!     GENERATION OF THE CONTINUOUS FREQUENCY GRID

      CALL FGRID(NFDIM,NF,XLAMBDA,FWEIGHT,AKEY,NOM,SYMBOL,NATOM,N,NCHARG,ELEVEL,EION,EINST)

      CALL GEOMESH(RADIUS, ENTOT, T, FAL, P, Z, RSTAR, AMU, ATMEAN, ND, NP)

      allocate(XJC(ND))
      allocate(XJCARR(ND, NF))
      allocate(XJL(ND, LASTIND))
      allocate(EDDI(3, ND))
      allocate(EDDARR(3, ND, NF))
      allocate(TAUROSS(ND))
      allocate(RNE(ND))
      allocate(VELO(ND))
      allocate(GRADI(ND))
      allocate(EMFLUX(NF))
      allocate(ENLTE(N))
      allocate(POPNUM(ND, N))
      allocate(HTOT(ND))
      allocate(GTOT(ND))
      allocate(XTOT(ND))
      allocate(ETOT(ND))
      allocate(POP1(ND, N))
      allocate(POP2(ND, N))
      allocate(POP3(ND, N))

      CALL mol_ab(ABXYZn, ABXYZ, SYMBOL, ENTOT, T, ND)

      IF(.NOT. ALLOCATED(ABXYZ_small))  allocate(ABXYZ_small(NATOM))
      IF(.NOT. ALLOCATED(ABXYZn_small)) allocate(ABXYZn_small(NATOM, ND))

      ABXYZ_small(1:NATOM)=ABXYZ(1:NATOM)

      do i = 1, ND

         do j = 1, NATOM

            ABXYZn_small(j,i) = ABXYZn(j,i)

         enddo

      enddo

!     RINAT TAGIROV:
!     Calculation or read-out of the velocity field.
!     In case of the calculation
!     the expansion law is regulated with the boundary
!     condition parameters VMIN and VFINAL
!     which are set in the CARDS file.
!     Using these parameters INITVEL subroutine
!     generates VPAR1 and VPAR2 parameters
!     entering the function WRVEL below which calculates the ultimate velocities.

!     In case of the read-out the velocity is read from the VEL_FIELD_FILE set in FILE_OPERATIONS.FOR
!     and if there is a local extremum in the input data...
!     (local maximum in case of negative velocities, local minimum in case of positive velocities;
!     at least it was like that in the 3D simulation data I worked with and => EXTRAP_VEL_FIELD subroutine
!     below works only for this two cases, however the rest of the code can't handle the negative case
!     because apparently the frequency initial condition in it is not suited for contracting flows,
!     see Mihalas, Kunasz & Hummer, ApJ, 202:465-489, 1975, Appendix B for details about various velocity laws and the
!     type of radiative transfer scheme implemented in the code)
!     ...the law gets extrapolated from the point of the extremum
!     to the innermost point yelding therefore a monotonically increasing/decreasing function.
!     The height grid in the VEL_FIELD_FILE has to be the same as in the atmosphere model file ATM_MOD.
!     The TABLE string in CARDS file was used before to
!     control the calculation/read-out option but is obsolete now (it is still in the CARDS file though).
!     The logical variable VEL_FIELD_FROM_FILE is declared in comblock.for and set in hminus.for.
!     The units of velocity in the VEL_FIELD_FILE are km/s, height is in km.
!     The first column is height, the second is velocity.

      IF (VEL_FIELD_FROM_FILE) THEN

         ALLOCATE(VELO_NE(ND))
         ALLOCATE(VELO_E(ND))

         FILE_UNIT = 1832

         OPEN(UNIT = FILE_UNIT, FILE = VEL_FIELD_FILE, ACTION = 'READ')

         DO I = 1, ND; READ(FILE_UNIT, *) H, VELO_NE(I); ENDDO

         CLOSE(FILE_UNIT)

!        Extrapolation of the velocity law. If there is no local minimum/maximum then VELO_E = VELO_NE.
         VELO_E = EXTRAP_VEL_FIELD(VELO_NE(1 : ND), ND)

         VELO(1 : ND) = VELO_E(1 : ND)

      ELSE

         DO L = 1, ND; VELO(L) = WRVEL(RADIUS(L)); ENDDO

      ENDIF

      CALL GRADIFF(ND,VELO,GRADI,RADIUS)
 
C***  STAPEL: NUMBER OF FREE ELECTRONS PER ATOM
C***  S T A R T   A P P R O X I M A T I O N

      STAPEL = 0.0d0

      DO NA = 1, NATOM; STAPEL = STAPEL + ABXYZ(NA) * (STAGE(NA) - 1.); ENDDO

      RNE(1 : ND) = STAPEL

C***  Read Line-blanketing table
      IF (LBLANK.NE.0) LBLANK=-2
      CALL REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)
      IF (ABS(LBLANK).EQ.2) 
     $CALL PRIBLA (LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,SCAFAC,ABSFAC)
 
C***  TEMPERATURE STRATIFICATION AND INITIAL POPNUMBERS (LTE)

      CALL GREYM(ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR,
     $           ALPHA,SEXPO,AGAUNT,POPNUM,TAUROSS,R23,TTABLE,
     $           LBKG,XLBKG1,XLBKG2,N,
     $           LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,KODAT,
     $           NOM,NFIRST,NLAST,NATOM,WAVARR,SIGARR)

      CALL PRIMOD(ND,RADIUS,RSTAR,ENTOT,T,VELO,GRADI,NP,MODHEAD,JOBNUM,TTABLE,TAUROSS,R23)

!=================================================================
! CONSTANT ELECTRON CONCENTRATION RUN AND/OR LTE RUN

      IF (LTE_RUN .OR. CONST_ELEC) THEN
!      IF (LTE_RUN) THEN

         ALLOCATE(ELEC_CONC(ND))
         ALLOCATE(HEAVY_ELEM_CONC(ND))

         ELEC_CONC =       READ_ATM_MOD(fal_mod_file, '3')
         HEAVY_ELEM_CONC = READ_ATM_MOD(fal_mod_file, '4')

         RNE(1 : ND) = ELEC_CONC(1 : ND) / HEAVY_ELEM_CONC(1 : ND)

         DEALLOCATE(ELEC_CONC)
         DEALLOCATE(HEAVY_ELEM_CONC)

      ENDIF

!=================================================================
 
!     MODFILE: MASS STORAGE FILE IN NAME-INDEX MODE
      TOTOUT =       0.0d0
      TOTIN =        0.0d0

      POP1(:, :) =   0.0d0
      POP2(:, :) =   0.0d0
      POP3(:, :) =   0.0d0

      HTOT(1 : ND) =   0.0d0
      GTOT(1 : ND) =   0.0d0
      XTOT(1 : ND) =   0.0d0
      ETOT(1 : ND) =   0.0d0

      EMFLUX(1 : NF) = 0.0d0

      if (allocated(wcharm)) deallocate(wcharm)

      allocate(wcharm(ND, NF))

      WCHARM(1 : ND, 1 : NF) = 0.0d0

      IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'UNKNOWN')

      JOBNUM = 0

      CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM)

      CLOSE(ifl)

      ifl = 3; open(ifl, file = 'POPNUM', status = 'unknown')

      call writpop(ifl, T, popnum, pop1, pop2, pop3, rne, n, nd, modhead, jobnum)

      close(ifl)

!     START APPROXIMATION FOR THE RADIATION FIELD
!     JSTART writes the files RADIOC and RADIOL

      EDDI(1 : 3, 1 : ND) = 0.0d0

      CALL JSTART(NF,XLAMBDA(1 : NF),ND,T,XJC,XJL,
     $            HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $            NCHARG,ELEVEL,EDDI,WCHARM,NOM,N,EINST,
     $            MODHEAD,JOBNUM,TEFF)

      write(*,  *) 'WRSTART - ', fstring, ' run time: ', TOC(timer)

      open(78, file = 'MODHIST', status = 'unknown')
   
      write(78, *) 'WRSTART - ', fstring, ' run time: ', TOC(timer)
   
      close(78)

      RETURN
 
      END SUBROUTINE


      FUNCTION EXTRAP_VEL_FIELD(VEL, ND) RESULT(VEL_E)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::               ND

      REAL*8, DIMENSION(ND), INTENT(IN) :: VEL
      REAL*8, DIMENSION(ND) ::             VEL_E

      INTEGER ::                           I, ML, NI, EXTRLOC

      REAL*8 ::                            A, B

      VEL_E(1 : ND) = VEL(1 : ND)

      IF (VEL(1) .GT. 0) EXTRLOC = MINLOC(VEL, 1)
      IF (VEL(1) .LT. 0) EXTRLOC = MAXLOC(VEL, 1)

!     If there is no local minimum/maximum then no extrapolation is needed
      IF (EXTRLOC .EQ. 1) THEN

         WRITE(*, '(/,A,1x,A,1x,A,/)') 'NO LOCAL MINIMUM/MAXIMUM WAS FOUND IN', VEL_FIELD_FILE,
     $                                 '=> NO EXTRAPOLATION HAS BEEN PERFORMED.
     $                                  PROCESSING THE VELOCITY FIELD AS IT IS.'

         RETURN

      ENDIF

      ML = EXTRLOC; NI = ML - 1

      A = (VEL(ML) - VEL(NI)) / (HEIGHT(ML) - HEIGHT(NI))

      B = VEL(ML) - A * HEIGHT(ML)

      DO I = ML + 1, ND, 1; VEL_E(I) = A * HEIGHT(I) + B; ENDDO

      RETURN

      END FUNCTION EXTRAP_VEL_FIELD

      END MODULE
