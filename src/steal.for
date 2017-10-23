      MODULE MOD_STEAL

      CONTAINS

      SUBROUTINE STEAL(job)

!     STATISTICAL EQUATIONS WITH APPROXIMATE LAMBDA-OPERATORS

      use MOD_LINPOP
      use MOD_PLOTPOP
      use MOD_PLOTT
      use MOD_POPZERO
      use MOD_PRIBLA
      use MOD_PRICORR
      use MOD_PRIDAT
      use MOD_PRIEXPO
      use MOD_PRIH
      use MOD_PRIHIST
      use MOD_PRIPOP
      use MOD_PRITAU
      use MOD_READMOD
      use MOD_READPOP
      use MOD_READRAD
      use MOD_REDCOR
      use MOD_STHIST
      use MOD_WRITPOP
      use MOD_PLOTFLU
      use MOD_PRIFLUX
      use MOD_REBLANK
      use MOD_WRITMOD

      use mod_decode
      use file_operations
      use common_block
      use vardatom_full
      use vardatom_nlte
      use varhminus
      use varsteal

      IMPLICIT REAL*8(A - H, O - Z)

CMH  CHANGES BY MARGIT HABERREITER
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
 
      COMMON /VELPAR/  VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE

      integer   NGAMR(10),NGAML(10)
      real*8    AGAMR(10),AGAML(10)

      logical   line(lastind_nlte)

      logical   TPLOT

      character MODHEAD*104,MODHOLD*104, CARD*80, LCARD*120

      character*7 JOB

      integer,external :: time

      integer :: tstart

      logical PROLIB

      REAL*8 :: CORMAX

      real*8 :: steal_start,  steal_finish

      call cpu_time(steal_start); call system("echo -n $(date +%s) >> wall_time.steal")

      WRITE(*, *), 'entered steal: job = '//JOB

      IF(LBKG) PRINT*, 'STEAL: LINE BLANKETING = TRUE'

      tstart = time()

      if (allocated(levelpl)) deallocate(levelpl); allocate(levelpl(N))

C***  DECODING INPUT DATA ******************************************
      CALL DECSTE(LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,IFRRA,ITORA,LSEXPO,
     $            IFLUX,IDAT,LEVELPL,N,IPLOTF,NEWWRC,
     $            NGAMR,NGAML,AGAMR,AGAML,LINE,lastind_nlte,TPLOT,
     $            Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB,LBLANK)

C***  READING OF THE MODEL FILE ----------------------------------------
      IFL = 3; open(IFL, file='MODFILE', STATUS='OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      if (allocated(itne))    deallocate(itne);    allocate(itne(ND))
      if (allocated(iwarn))   deallocate(iwarn);   allocate(iwarn(ND))

      if (allocated(depart))  deallocate(depart);  allocate(depart(ND, N))

      IFL = 3; open(IFL, file='POPNUM', STATUS='OLD')

c***  pop1 is dummy read because it will be overwritten below
      call readpop(ifl, T, popnum, pop2, pop3, pop1, rne, n, nd, modhead, jobnum)

      close(IFL)

c***  read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)	

      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             ncharg_nlte,EDDARR,EDDI,nom_nlte,WCHARM,
     $             N_nlte,lastind_nlte,einst_nlte,MODHEAD,JOBNUM)

!      stop 'steal stop'

c***  advance job-number counter
      JOBNUM = JOBNUM + 1

      CALL change(popnum, pop1, nd*n)

      CALL REBLANK(LBLANK, NF, XLAMBDA, ND, ENTOT, RNE, SCAFAC, ABSFAC)

      if (lblank .lt. 0) then
c***     the new blanketing table needs to be written to the model file
         IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'UNKNOWN')

         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $                XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)

         CLOSE(ifl)

      endif

      GAMMAR = 0.0D0
      GAMMAL = 0.0D0
 
      IF (JOBNUM .LE. 1) THEN

         CALL POPZERO(T,
     $                RNE,
     $                POPNUM,
     $                DEPART,
     $                ENTOT,
     $                ITNE,
     $                N,
     $                WEIGHT,
     $                NCHARG,
     $                EION,
     $                ELEVEL,
     $                EINST,
     $                LEVEL,
     $                FWEIGHT,
     $                XJCARR,
     $                NF,
     $                XJL,
     $                IFRRA,
     $                ITORA,
     $                AGAUNT,
     $                MODHEAD,
     $                MODHOLD,
     $                JOBNUM,
     $                ND,
     $                LSRAT,
     $                ALTESUM,
     $                COCO,
     $                KEYCOL,
     $                NOM,
     $                NATOM,
     $                KODAT,
     $                NFIRST,
     $                NLAST)

      ELSE

!     RATE EQUATION WITH APPROXIMATE RADIATION TRANSFER
!     IN THIS BRANCH, PRIRAT MAY ONLY SHOW THE NETTO RATES
!     CALCULATION OF NEW POPULATION NUMBERS, EL. DENSITY AND DEPARTURE COEFF.

         CALL LINPOP(T,
     $               RNE,
     $               ENTOT,
     $               ITNE,
     $               POPNUM,
     $               DEPART,
     $               POP1,
     $               N_nlte,
     $               weight_nlte,
     $               ncharg_nlte,
     $               eion_nlte,
     $               elevel_nlte,
     $               einst_nlte,
     $               level_nlte,
     $               XLAMBDA,
     $               FWEIGHT(1 : NF),
     $               XJCARR,
     $               NF,
     $               NFDIM,
     $               XJL,
     $               WCHARM,
     $               EPSILON,
     $               MODHEAD,
     $               JOBNUM,
     $               IFRRA,
     $               ITORA,
     $               RADIUS,
     $               RSTAR,
     $               IWARN,
     $               MAINPRO,
     $               MAINLEV,
     $               VDOP,
     $               indnup_nlte,
     $               indlow_nlte,
     $               lastind_nlte,
     $               ND,
     $               LSRAT,
     $               alpha_nlte,
     $               sexpo_nlte,
     $               agaunt_nlte,
     $               coco_nlte,
     $               keycol_nlte,
     $               altesum_nlte,
     $               nom_nlte,
     $               natom_nlte,
     $               kodat_nlte,
     $               levatnum_nlte,
     $               nfirst_nlte,
     $               nlast_nlte,
     $               WAVARR,
     $               SIGARR,
     $               JOBMAX,
     $               N,
     $               weight,
     $               ncharg,
     $               eion,
     $               elevel,
     $               einst,
     $               level,
     $               alpha,
     $               sexpo,
     $               agaunt,
     $               nom,
     $               natom,
     $               nfirst,
     $               nlast)

      ENDIF
 
C***  REDUCED CORRECTIONS, IF OPTION IS SET
      IF (REDUCE .NE. 1.0D0 .AND. JOBNUM .GT. 1) CALL REDCOR(POPNUM,POP1,ND,N,RNE,NCHARG,REDUCE)
 
      IF (LSPOP.GT.0) CALL PRIPOP(LSPOP,WEIGHT,NCHARG,NOM,ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)
 
      IF (JOBNUM .GT. 1) CALL PRICORR(POPNUM,POP1,LEVEL,N,ND,MODHEAD,LSPOP,CORMAX,
     $                                NCHARG,RNE,JOBNUM,REDUCE,GAMMAL,GAMMAR,
     $                                T,ELEMENT,NATOM,NFIRST,NLAST)

      IF (LSEXPO .GT. 0.AND.JOBNUM.GT.3) CALL PRIEXPO(POPNUM,POP1,POP2,LEVEL,N,ND,MODHEAD,JOBNUM,LSEXPO)
 
C***  UPDATING THE MODEL FILE
      ifl = 2; OPEN(ifl, FILE = 'POPNUM', STATUS = 'unknown')

      call writpop(ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum); CLOSE(ifl)

c***  if a new LB table is read then store a new model file
      if (lblank.lt.0) then
      IFL=3
      open (IFL,file='MODFILE',STATUS='UNKNOWN')
      CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM)
      CLOSE (ifl)
      endif

C***  FIND LASTWRC = JOBNUMBER OF LAST WRCONT JOB
      LASTWRC=0
      open (79,file='MODHIST',status='old')
  8   read (79,'(A80)',end=11) card
      do while (card(7:12).ne.'WRCONT')
	   read (79,'(A80)',end=11,err=555) card
      enddo
      read (card,'(1X,I3)') LASTWRC
	goto 8
 11   continue
C***  UPDATING THE MODEL HISTORY
      CALL STHIST(LCARD,GAMMAL,GAMMAR,MODHEAD,JOBNUM,CORMAX,REDUCE,MODHOLD,time()-tstart)
      write (79,'(A120)') LCARD
      close (79)

      IF (JOBNUM .GE. JOBMAX ) THEN
         WRITE (6,*) 'MAX. NUMBER OF JOBS EXCEEDED'
         PRINT *,' MAX. NUMBER OF JOBS EXCEEDED'
         GOTO 20
         ENDIF
 
      IF (JOBNUM .EQ. 1) GOTO 15

      IF (CORMAX .LT. EPSILON) THEN

         PRINT *,' REPEAT CYCLE IS CONVERGED'
         JOBDIFF=JOBNUM-LASTWRC
         IF (JOBDIFF .LT. 0) JOBDIFF=JOBDIFF+100
         IF (JOBDIFF .LE. 3) THEN
            PRINT *,' -------  MODEL FINALLY CONVERGED !  ---------'
            IHIST=1
            IFLUX=1
            IDAT=1
            IF (LBLANK.EQ.1) LBLANK=2
            LSPOP=1
            LPRIH=1

            CALL PRIPOP(LSPOP,WEIGHT,NCHARG,NOM,ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)

            CALL PRITAU(MODHEAD,JOBNUM,RSTAR,ND,RADIUS,RNE,ENTOT,T,
     $                  POPNUM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,
     $                  EION,EINST,ALPHA,SEXPO,AGAUNT,NOM,XLAMBDA,
     $                  FWEIGHT,TAUROSS,WAVARR,SIGARR,NF,NFDIM)

            CALL PRIH(LPRIH,ND,RADIUS,HTOT,TEFFE,T,TAUROSS,JOBNUM,MODHEAD)
            GOTO 20
            ELSE
            IF (LSEXPO .NE. 1)
     $       CALL PRIEXPO (POPNUM,POP1,POP2,LEVEL,N,ND,MODHEAD,JOBNUM,1)
              GOTO 15

          ENDIF

      ELSE ! not converged

        WRITE(*, '(A)') 'STEAL: LAMBDA ITERATION CYCLE HAS NOT CONVERGED YET'

        IF (JOBNUM - LASTWRC .GT. 3 * (NEWWRC - 1)) THEN

          WRITE(*, '(A)'), 'STEAL: EDDINGTON FACTORS ARE TOO OLD'; GOTO 15

        ELSE

          IF (JOBNUM .LE. 25) THEN

              WRITE (*, '(A)') 'STEAL: FIRST 6 ITTERATIONS => NEW EDDINGTON FACTORS WILL BE CALCULATED'; GOTO 15

          ENDIF

          WRITE(*, '(A)') 'STEAL: NEW REPEAT JOB TO BE ROUTED'; REWIND 99; WRITE(99, '(A6)') 'repeat'; JOB = 'repeat'; GOTO 30

        ENDIF

      ENDIF ! converged / not converged
 
!     BRANCH FOR ROUTING WRCONT OR EXTRAP JOB
   15 IF ((JOBNUM .LE. 30.) .OR. (CORMAX .LT. EPSILON)) THEN

         WRITE(*, '(A)'), 'STEAL: WRCONT TO BE ROUTED'; REWIND 99; WRITE(99, '(A6)') 'wrcont'; JOB = 'wrcont'

      ELSE

         WRITE(*, '(A)'), 'STEAL: EXTRAP IS NEXT'; REWIND 99; WRITE(99, *) 'extrap'; JOB = 'extrap'

      ENDIF
 
      GOTO 30
 
C***  BRANCH FOR NO SUBSEQUENT JOB
   20 WRITE(*, '(A)'), 'STEAL: NO SUBSEQUENT JOB TO BE ROUTED'

      REWIND 99; WRITE(99, '(a4)') 'exit'; JOB = 'exit'

C***  PRINTOUT OF MODEL HISTORY (IF REQUESTED)
      IF (IHIST.EQ.1) CALL PRIHIST(MODHEAD,JOBNUM)
C***  PRINTOUT OF EMERGENT CONT. FLUX (IF REQUESTED)
      IF (IFLUX .EQ. 1) THEN
         IF (TOTOUT .EQ. 0.d0 ) THEN
            PRINT 7
            ELSE
            CALL PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,JOBNUM,
     $                    FWEIGHT,MODHEAD,AKEY )
            ENDIF
         ENDIF
C***  PRINTOUT OF ATOMIC DATA (IF REQUESTED)
      IF (IDAT.EQ.1)
     $CALL PRIDAT(N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,EINST,
     $            KODAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $            NATOM,ELEMENT,NOM,ABXYZ,ATMASS)

      IF (ABS(LBLANK).EQ.2) 
     $CALL PRIBLA (LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,
     $	                 SCAFAC,ABSFAC)
C***  DIRECT TRANSFER OF EMERGENT CONT. FLUX PLOT (IF REQUESTED)
      IF (IPLOTF .EQ. 1) THEN
         IF (TOTOUT .EQ. 6HUNDEF. ) THEN
            PRINT 7
    7       FORMAT (//' INVALID PRINT OR PLOT OPTION - ',
     $      'EMERGENT CONT. FLUX NOT YET CALCULATED ',//)
            ELSE
C***  Y0 IST DER DISTANZMODULUS
C***  -20.00471 IST LOG10(PI*C/(10PARCEC)**2) IN A/CM2)
      Y0=-Y0/2.5+2.*LOG10(RSTAR)-20.00471
            CALL PLOTFLU (NF,XLAMBDA,EMFLUX,MODHEAD,JOBNUM,
     $           Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB)
            ENDIF
         ENDIF
C***  DIRECT TRANSFER OF POPNUMBER PLOT (IF REQUESTED)
      IF (LEVELPL(1) .NE. 0)
     $      CALL PLOTPOP (LEVELPL,N,ND,LEVEL,ENTOT,POPNUM,MODHEAD,
     $          JOBNUM )
C***  DIRECT TRANSFER OF TEMPERATURE STRATIFICATION PLOT (IF REQUESTED)
      IF (TPLOT) CALL PLOTT (ND,RADIUS,T,MODHEAD,JOBNUM)
 
C***  PROGRAM STOP
   30 CONTINUE

      CLOSE(99)

      call system("echo ' '$(date +%s) >> wall_time.steal")

      call cpu_time(steal_finish)

      call open_to_append(302, 'cpu_time.steal'); write(302, '(F6.3)') steal_finish - steal_start; close(302)

      return

555   continue
      write (6,*) ' error during history file read - search for wrcont'
      stop 'error ft7'

666   continue
      write (6,*) ' error during history file read - forward to EOF'
      stop 'error ft7'

      end subroutine

      subroutine change(array1, array2, N)

      implicit none

      integer, intent(in) ::               N

      real*8, dimension(N), intent(in) ::  array1

      real*8, dimension(N), intent(out) :: array2

      array2(1 : N) = array1(1 : N)

      return

      end subroutine

      end module
