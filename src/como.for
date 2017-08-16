      module MOD_COMO

      contains

      SUBROUTINE COMO

      use MOD_CALCLAMBDAS
      use MOD_COOP_M
      use MOD_DATOM_M
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
      use PARAMS_ARRAY ! NFDIM is known from here (see params_array.for)
      use ABUNDANCES
      USE MOD_BFCROSS

      use common_block
      use vardatom
      use varhminus
      use varsteal
      use file_operations

      IMPLICIT NONE

      integer, parameter:: IPDIM = 25, NBDIM=99
      
      integer IFL, IPMAX
      integer JOBNUM
      integer K, KEYCON
      integer LASTIND,LASTK,LBLANK,LBLAON,LSINT,LSOPA

      integer MAXVAL,MINVAL
      integer N,NATOM,NBINW,NBMAX,NCON,ND

      integer NF,NFCDIM,NP
      real*8 A, ABSEVT, ABSFAC, ALMAX, ALMIN
      real*8 B
      real*8 C, CARD, CRATE, DM
      real*8 HNU
      real*8 RSTAR
      real*8 SCAEVT,SCAFAC,SCAGRI
      real*8 TEFF,TOTIN
      real*8 TOTOUT,U,VDOP,VJL,VL
      real*8 W
      real*8,allocatable:: DUMMY1(:)
      character*8,allocatable :: CDUMMY1(:)
      integer tdiff,tstart,tend
      integer,external :: time

      integer :: iii

!     CONTINUOUS RADIATION TRANSFER (MOMENT EQUATIONS) WITH GIVEN EDDI-FACTORS
!     FORMAL SOLUTION FROM GIVEN POP NUMBERS

      PARAMETER (NFCDIM = 40)
 
      COMMON // A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM)
     $ ,U(NDDIM,NPDIM),VL(NPDIM),HNU(NDDIM),VJL(NPDIM,NDDIM)
     $ ,KONOPT(NDIM),KEYCON(NFDIM)

      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)
      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2 
      CHARACTER MODHEAD*104, LCARD*100

      CHARACTER*10 MAINPRO(NDDIM),MAINLEV(NDDIM)
      character*80 konopt

CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF

      INTEGER XLBKG1,XLBKG2
      LOGICAL LBKG

      real*8 :: start, finish

      call cpu_time(start); call system("echo -n $(date +%s) >> wall_time.como")

      print*, 'entering como...'

      tstart = time()

      CALL DATOM_M(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $             EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $             INDNUP,INDLOW,LASTIND,NATOM,
     $             ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $             NLAST,WAVARR,SIGARR,NFDIM)

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

!      print*, 'como shape print: ', shape(popnum), sum(popnum(64, :)) + rne(64)

!      do iii = 1, 12

!         print*, 'como popnum print: ', iii, popnum(64, iii)

!      enddo

!      stop

      if (allocated(wcharm)) deallocate(wcharm)

      allocate(WCHARM(ND, NF))

      WCHARM(1 : ND, 1 : NF) = 0.0d0

!     read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             NCHARG,EDDARR,EDDI,NOM,WCHARM,N,lastind,
     $             EINST,MODHEAD,JOBNUM)

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

!         print*, 'keycon(k) = ', KEYCON(K) .EQ. 4HETLA

!     ONLY FREQUENCIES WHICH ARE NOT YET TREATED BY ETLA
!     IF "PRINT OPA" OPTION IS GIVEN, THE JUMP IS AFTER "CALL PRIOPA"
         IF (KEYCON(K).EQ.4HETLA .AND. LSOPA.LE.0 ) then

             print*, 'check 1 is true'

             GOTO 6

         endif
 
         CALL COOP_M(XLAMBDA(K),ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $               OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $               N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $               DUMMY1,DUMMY1,CDUMMY1,K,SIGMAKI,
     $               WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF),
     $               LBKG,XLBKG1,XLBKG2,NF)

         IF (LSOPA .GT. 0) THEN

             CALL PRIOPA(XLAMBDA(K),K,ND,LSOPA,RADIUS,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)

         ENDIF

         IF (KEYCON(K) .EQ. 4HETLA) then

            print*, 'check 2 is true'

            GOTO 6

         endif
 
!        now extract XJC and EDDI for the frequency K
         call extrxjc(XJCARR,XJC,EDDARR,EDDI,nd,nf,K)

         CALL MOMO(OPA,ETA,THOMSON,EDDI,RADIUS,XJC,A,B,C,W,ND)

         IF (LSINT.GT.0) CALL PRIMINT(XJCARR,ND,XLAMBDA,NF,K,LSINT,EDDI,JOBNUM,MODHEAD)

   6     CONTINUE

         WCHARM(1 : ND, K) = CALCLAMBDAS(OPA, RADIUS, EDDI, ND)

C***     UPDATING THE CONTINUOUS RADIATION FIELD ON THE MODEL FILE
c***     XJC and EDDI are stored for later write to file RADIOC

         call storxjc(XJCARR, XJC, EDDARR, EDDI, nd, nf, K)
          
!         IF (LTE_RUN) XJC_LTE(1 : ND, K) = XJCARR(1 : ND, K)
!
!         IF (LTE_RUN) CALL PRINT_LTE_CONT(XLAMBDA(K), K, XJC_LTE(1 : ND, K))

      ENDDO

!      stop

      call assert(.not.any(isnan(WCHARM)),'COMO: WCHARM is NaN')

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

!      do k = 1, NF

!         do iii = 1, ND

!            write(*, '(A,1x,i4,1x,e15.7,1x,i4,1x,e15.7)'), 'como wcharm: ', k, xlambda(k), iii, wcharm(iii, k)

!         enddo

!      enddo

!      stop 'como stop'
     
C***  ENDLOOP  ---------------------------------------------------------
 
C***  NOTE THAT THE POPNUMBERS ARE NOT UPDATED BY THIS PROGRAM  !!!!!
      call writradc(xjcarr,xjc,eddarr,eddi,emflux,totin,totout,
     $              HTOT,GTOT,XTOT,ETOT,WCHARM,nd,nf,MODHEAD,JOBNUM)
 
C***  UPDATING THE MODEL HISTORY
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

      call system("echo ' '$(date +%s) >> wall_time.como"); call cpu_time(finish)

      call open_to_append(2313, 'cpu_time.como'); write(2313, '(F6.3)') finish - start; close(2313)

      return

      end subroutine

      end module
