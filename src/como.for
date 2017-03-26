      module MOD_COMO

      contains

      SUBROUTINE COMO

      use MOD_COOP
      use MOD_DATOM
      use MOD_DECOMO
      use MOD_DERIV
      use MOD_ERROR
      use MOD_extrxjc
      use MOD_MOMO
      use MOD_PRIBLA
      use MOD_PRIMINT
      use MOD_PRIOPA
      use MOD_READMOD
      use MOD_READPOP
      use MOD_READRAD
      use MOD_REBLANK
      use MOD_STORXJC
      use MOD_WRITMOD
      use MOD_WRITRADC
      use UTILS
      use ABUNDANCES
      USE MOD_BFCROSS

      use common_block
      use vardatom
      use varhminus
      use varsteal
      use phys
      use local_operator
      use file_operations

      implicit none

      integer IFL
      integer JOBNUM
      integer K, KEYCON
      integer LASTIND,LASTK,LBLANK,LSINT,LSOPA

      integer MAXVAL,MINVAL
      integer N,NATOM,NCON,ND

      integer NF, NP
      real*8 CARD, CRATE, DM
      real*8 RSTAR
      real*8 TEFF,TOTIN
      real*8 TOTOUT,VDOP
      real*8,allocatable:: DUMMY1(:)
      character*8,allocatable :: CDUMMY1(:)
      integer tdiff,tstart,tend
      integer,external :: time

      integer :: l

!     CONTINUOUS RADIATION TRANSFER (MOMENT EQUATIONS) WITH GIVEN EDDI-FACTORS
!     FORMAL SOLUTION FROM GIVEN POP NUMBERS

      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2 
      CHARACTER MODHEAD*104, LCARD*100

      character*80 konopt(2)

CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF

      INTEGER XLBKG1,XLBKG2
      LOGICAL LBKG

      tstart = time()

      CALL DATOM(datom_lte,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $           EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $           INDNUP,INDLOW,LASTIND,NATOM,
     $           ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $           NLAST,WAVARR,SIGARR,NFDIM)

      IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'OLD')

      CALL DECOMO(LSOPA, LSINT, KONOPT, NCON, LBLANK)

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      IFL = 3; open(IFL, file = 'POPNUM', STATUS = 'OLD')

      call readpop(ifl, T, popnum, pop1, pop2, pop3, rne, n, nd, modhead, jobnum)

      close(ifl)

      if (allocated(wcharm))    deallocate(wcharm)

      if (allocated(tau_cont))  deallocate(tau_cont)

      if (allocated(damp_cont)) deallocate(damp_cont)

      allocate(WCHARM(ND, NF))

      allocate(tau_cont(ND, NF))

      allocate(damp_cont(ND, NF))

      damp_cont(1 : ND, 1 : NF) = .false.

      WCHARM(1 : ND, 1 : NF) = 0.0d0

!     read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             ncharg_nlte,EDDARR,EDDI,nom_nlte,WCHARM,N_nlte,lastind_nlte,
     $             einst_nlte,MODHEAD,JOBNUM)

      JOBNUM = JOBNUM + 1

!     the blanketing table is read by routine READMOD
!     if lblank.gt.0 then read a new table from the file LIBLANK
      CALL REBLANK(LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

      if (lblank .lt. 0) then
!     the new blanketing table needs to be written to the model file
         IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'UNKNOWN')

         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,
     $                XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)

         close(ifl)

      endif

      IF (abs(LBLANK) .EQ. 2) CALL PRIBLA(LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,SCAFAC,ABSFAC)
 
!     PRECALCULATION OF THE BOUND-FREE CROSS SECTIONS SIGMAKI
      CALL BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,
     $             XLAMBDA(1 : NF),ALPHA,SEXPO,AGAUNT,NOM,
     $             WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF))

!     ETLA TREATMENT OF CONTINUA (OPTIONALLY)
      IF (NCON .EQ. 0) GOTO 4
 
    4 CONTINUE
 
!     LOOP OVER ALL FREQUENCY POINTS
!     SOLUTION OF THE MOMENT EQUATION AT EACH FREQUENCY POINT (FORMAL SOLUTION)
!     WRITE OUTPUT ETA,OPA TO FILE

      DO K = 1, NF

         lastk = K

!     ONLY FREQUENCIES WHICH ARE NOT YET TREATED BY ETLA
!     IF "PRINT OPA" OPTION IS GIVEN, THE JUMP IS AFTER "CALL PRIOPA"
         IF (KEYCON(K).EQ.4HETLA .AND. LSOPA.LE.0 ) GOTO 6
 
         CALL COOP(XLAMBDA(K),ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $             OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $             N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $             DUMMY1,DUMMY1,CDUMMY1,K,SIGMAKI,
     $             WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF),
     $             LBKG,XLBKG1,XLBKG2,NF)

!         print*, 'como: K = ', K, ' out of ', NF

         IF (LSOPA .GT. 0) THEN

             CALL PRIOPA(XLAMBDA(K),K,ND,LSOPA,RADIUS,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)

         ENDIF

         IF (KEYCON(K) .EQ. 4HETLA) GOTO 6
 
!        now extract XJC and EDDI for the frequency K
         call extrxjc(XJCARR,XJC,EDDARR,EDDI,nd,nf,K)

         CALL MOMO(OPA,ETA,THOMSON,EDDI,RADIUS,XJC,ND)

         IF (LSINT.GT.0) CALL PRIMINT(XJCARR,ND,XLAMBDA,NF,K,LSINT,EDDI,JOBNUM,MODHEAD)

   6     CONTINUE

         WCHARM(1 : ND, K) = cont_loc_oper(OPA, RADIUS, EDDI, ND)

         tau_cont(1 : ND, K) = opt_dep(opa(1 : ND) / RSTAR, 1.0D+5 * height(1 : ND), ND)

!        UPDATING THE CONTINUOUS RADIATION FIELD ON THE MODEL FILE
!        XJC and EDDI are stored for later write to file RADIOC

         call storxjc(XJCARR, XJC, EDDARR, EDDI, nd, nf, K)
          
!         IF (LTE_RUN) XJC_LTE(1 : ND, K) = XJCARR(1 : ND, K)
!
!         IF (LTE_RUN) CALL PRINT_LTE_CONT(XLAMBDA(K), K, XJC_LTE(1 : ND, K))

      ENDDO

      call assert(.not.any(isnan(WCHARM)),'COMO: WCHARM is NaN')

!     perform the acceleration damping in case the density is too high at the outer edge of the atmosphere
      if (damp_acc) call acc_damp(ND, NF, tau_cont, WCHARM, damp_cont)

      !Sanity Check ----------------------------------------------------
      IF(maxval(WCHARM) >= 1d0-1d-20) THEN

         print '("como: WARN:max(WCHARM)=",e10.4,", RESET")',maxval(WCHARM)
         print*, 'location: ', maxloc(wcharm)
!         where(WCHARM >= 1d0-1d20 ) WCHARM = 1d0-1d20

      ENDIF

      IF(minval(WCHARM) < 1d-35) THEN

         print '("como: WARN:min(WCHARM)=",e10.4,", RESET")',minval(WCHARM)
         print*, 'location: ', minloc(wcharm)
!         where(WCHARM <  1d-35 ) WCHARM = 1d-35

      ENDIF

      call print_clo(ND, NF, xlambda, T, tau_cont, wcharm, damp_cont)

!     ENDLOOP  ---------------------------------------------------------
 
!     NOTE THAT THE POPNUMBERS ARE NOT UPDATED BY THIS PROGRAM  !!!!!
      call writradc(xjcarr,xjc,eddarr,eddi,emflux,totin,totout,
     $              HTOT,GTOT,XTOT,ETOT,WCHARM,nd,nf,MODHEAD,JOBNUM)
 
!     UPDATING THE MODEL HISTORY
      open (7,file='MODHIST',status='old')
  8   read (7,'(A80)',end=11) card
      goto 8
 11   continue
      tend=time()
      tdiff = tend-tstart
      write (LCARD,88) JOBNUM,'. COMO   1 -',LASTK,tdiff,' sec'
   88 FORMAT (1H/,I3,A,I3,' COMPLETE  -  run time all K: ',F10.2,A)
      write (7,'(A100)') LCARD
      close (7)   ! MODHIST

      RETURN

      END SUBROUTINE

      subroutine print_clo(nd, nf, wvl, T, tau, clo, damp)

      use common_block
      use file_operations

      implicit none

      integer, intent(in) :: nd, nf

      real*8, intent(in), dimension(nf) :: wvl

      real*8, intent(in), dimension(nd) :: T

      real*8, intent(in), dimension(nd, nf) :: tau, clo

      logical, intent(in), dimension(nd, nf) :: damp

      character(len = 1000) :: fmt_head, fmt_body

      integer :: file_unit

      character(len = 3) :: file_name

      integer :: l, k

      file_unit = 165; file_name = 'CLO'

      fmt_head = '(1x,A,4x,A,8x,A,8x,A,6x,A,8x,A,12x,A,14x,A,8x,A,/)'

      fmt_body = '(2(i5,2x),e15.7,2x,i3,2(3x,F9.2),2(2x,es15.7),4x,L)'
      
      if ((each_ali .and. lambda_iter .eq. 0) .or. .not. each_ali) then

          call rm_file(file_name, '-f')

          call open_to_append(file_unit, file_name)

          write(file_unit, fmt_head) 'lit', 'wid', 'wvl', 'hid', 'hei', 'tem', 'tau', 'clo', 'dam'

      else

          call open_to_append(file_unit, file_name)

      endif

      do k = 1, nf

         do l = 1, nd

            write(file_unit, fmt_body) lambda_iter, k, wvl(k), l, height(l), T(l), tau(l, k), clo(l, k), damp(l, k)

         enddo

      enddo

      close(file_unit)

      end subroutine print_clo

      end module
