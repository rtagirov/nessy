      MODULE MOD_STEAL

      CONTAINS

      SUBROUTINE STEAL(job)

!     STATISTICAL EQUATIONS WITH APPROXIMATE LAMBDA-OPERATORS

      use MOD_DATOM
      use MOD_DECSTE
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
      use ABUNDANCES

      use file_operations
      use common_block
      use vardatom
      use varhminus
      use varsteal

      IMPLICIT REAL*8(A - H, O - Z)

CMH  CHANGES BY MARGIT HABERREITER
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
 
      COMMON /VELPAR/  VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE

      COMMON /COMLBKG/ LBKG, XLBKG1, XLBKG2

      integer   NGAMR(10),NGAML(10)
      real*8    AGAMR(10),AGAML(10)

      integer   XLBKG1, XLBKG2

      logical   line(lastind_nlte)

      logical   TPLOT,NODM,LBKG

      character MODHEAD*104,MODHOLD*104, CARD*80, LCARD*120

      character*7 JOB

      integer,external :: time

      integer :: tstart

      logical PROLIB

      REAL*8 :: CORMAX

      WRITE(*, *), 'entering steal: job = '//JOB

      IF(LBKG) PRINT*, 'STEAL: LINE BLANKETING = TRUE'

      tstart = time()

C***  READING THE ATOMIC DATA FROM FILE DATOM

	CALL DATOM(datom_lte,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $         EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $         INDNUP,INDLOW,LASTIND,NATOM,
     $         ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $         NLAST,WAVARR,SIGARR,NFDIM)

      if (allocated(levelpl)) deallocate(levelpl); allocate(levelpl(N))
      if (allocated(nfedge))  deallocate(nfedge);  allocate(nfedge(N))

C***  DECODING INPUT DATA ******************************************
      CALL DECSTE(LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,IFRRA,ITORA,LSEXPO,
     $            IFLUX,IDAT,LEVELPL,N,IPLOTF,NEWWRC,
     $            NGAMR,NGAML,AGAMR,AGAML,DELTAC,LINE,lastind_nlte,TPLOT,
     $            Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB,LBLANK)

C***  READING OF THE MODEL FILE ----------------------------------------
      IFL = 3; open(IFL, file='MODFILE', STATUS='OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      if (allocated(opa))     deallocate(opa);     allocate(opa(ND))
      if (allocated(eta))     deallocate(eta);     allocate(eta(ND))

      if (allocated(thomson)) deallocate(thomson); allocate(thomson(ND))
      if (allocated(tauthom)) deallocate(tauthom); allocate(tauthom(ND))
      if (allocated(itne))    deallocate(itne);    allocate(itne(ND))
      if (allocated(iwarn))   deallocate(iwarn);   allocate(iwarn(ND))

      if (allocated(opac))    deallocate(opac);    allocate(opac(NF))
      if (allocated(etac))    deallocate(etac);    allocate(etac(NF))
      if (allocated(dopa))    deallocate(dopa);    allocate(dopa(NF))
      if (allocated(deta))    deallocate(deta);    allocate(deta(NF))
      if (allocated(expfac))  deallocate(expfac);  allocate(expfac(NF))

      if (allocated(sigmaki)) deallocate(sigmaki); allocate(sigmaki(NF, N))

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

         CALL POPZERO(T,RNE,POPNUM,DEPART,ENTOT,ITNE,N,ENLTE,
     $                WEIGHT,NCHARG,EION,ELEVEL,EINST,LEVEL,
     $                XLAMBDA,FWEIGHT,XJCARR,NF,XJL,IFRRA,ITORA,ALPHA,
     $                SEXPO,AGAUNT,MODHEAD,MODHOLD,JOBNUM,
     $                ND,LSRAT,SIGMAKI,ALTESUM,COCO,KEYCOL,NOM,NATOM,
     $                KODAT,NFIRST,NLAST,WAVARR,SIGARR)

      ELSE

!     RATE EQUATION WITH APPROXIMATE RADIATION TRANSFER
!     IN THIS BRANCH, PRIRAT MAY ONLY SHOW THE NETTO RATES
!     CALCULATION OF NEW POPULATION NUMBERS, EL. DENSITY AND DEPARTURE COEFF.

         CALL LINPOP(T,RNE,ENTOT,ITNE,POPNUM,DEPART,POP1,
     $               N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EINST,LEVEL,
     $               XLAMBDA,FWEIGHT(1 : NF),XJCARR,NF,XJL,WCHARM,
     $               EPSILON,NODM,DELTAC,MODHEAD,JOBNUM,IFRRA,ITORA,
     $               RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $               VELO,GRADI,VDOP,INDNUP,INDLOW,LASTIND,
     $               OPAC,DOPA,DETA,SIGMAKI,ND,LSRAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,
     $               ALTESUM,ETAC,NFEDGE,EXPFAC,NOM,NATOM,KODAT,NFIRST,
     $               NLAST,WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,JOBMAX)

      ENDIF
 
C***  REDUCED CORRECTIONS, IF OPTION IS SET
      IF (REDUCE .NE. 1.0D0 .AND. JOBNUM .GT. 1) CALL REDCOR(POPNUM,POP1,ND,N,RNE,NCHARG,REDUCE)
 
      IF (LSPOP.GT.0) CALL PRIPOP(LSPOP,WEIGHT,NCHARG,NOM,ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)
 
      IF (JOBNUM .GT. 1) CALL PRICORR(POPNUM,POP1,LEVEL,N,ND,MODHEAD,LSPOP,CORMAX,
     $                                NCHARG,RNE,JOBNUM,REDUCE,GAMMAL,GAMMAR,
     $                                T,DELTAC,ELEMENT,NATOM,NFIRST,NLAST)

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
      CALL STHIST(LCARD,GAMMAL,GAMMAR,DELTAC,MODHEAD,JOBNUM,CORMAX,REDUCE,MODHOLD,time()-tstart)
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
     $                  FWEIGHT,TAUTHOM,TAUROSS,WAVARR,SIGARR,NF)

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

      RETURN

555   continue
      write (6,*) ' error during history file read - search for wrcont'
      stop 'error ft7'

666   continue
      write (6,*) ' error during history file read - forward to EOF'
      stop 'error ft7'

      END SUBROUTINE

      SUBROUTINE CHANGE(array1, array2, N)
C***  THIS SUBROUTINE copies an array

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION array1(N), array2(N)
      DO 1 I = 1, N

         array2(i) = array1(i)

    1 CONTINUE

      RETURN

      END subroutine

      END MODULE
