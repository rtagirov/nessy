      MODULE MOD_STEAL

      CONTAINS

      SUBROUTINE STEAL (job)

!     STATISTICAL EQUATIONS WITH APPROXIMATE LAMBDA-OPERATORS

      use MOD_CHANGE
      use MOD_DATOM_M
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

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      use VARDATOM

      IMPLICIT REAL*8(A - H, O - Z)

C******************************************************************************
CMH  CHANGES BY MARGIT HABERREITER
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
C******************************************************************************
 
      COMMON // RADIUS(NDDIM),ENTOT(NDDIM),T(NDDIM)
     $ ,XJC(NDDIM),XJCARR(NDDIM,NFDIM),XJL(NDDIM,MAXIND)
     $ ,EDDI(3,NDDIM),EDDARR(3,NDDIM,NFDIM),TAUROSS(NDDIM)
     $ ,RNE(NDDIM),VELO(NDDIM),GRADI(NDDIM),AINCRIT(NDDIM)
     $ ,XLAMBDA(NFDIM),FWEIGHT(NFDIM),EMFLUX(NFDIM),AKEY(NFDIM)
!     $ ,WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
!     $ ,ELEVEL(NDIM),EION(NDIM)
!     $ ,EINST(NDIM,NDIM),ALPHA(NDIM),SEXPO(NDIM)
     $ ,ENLTE(NDIM)
!     $ ,COCO(NDIM,NDIM,4),ALTESUM(4,NDIM)
     $ ,P(NPDIM),Z(NDDIM,NPDIM),POPNUM(NDDIM,NDIM)
!     $ ,ATMASS(MAXATOM),STAGE(MAXATOM),AGAUNT(NDIM)
     $ ,HTOT(NDDIM),GTOT(NDDIM),XTOT(NDDIM),ETOT(NDDIM)
     $ ,POP1(NDDIM,NDIM),POP2(NDDIM,NDIM),POP3(NDDIM,NDIM)
     $ ,EN(NDIMP2),V1(NDIMP2),V2(NDIMP2),RATCO(NDIMP2,NDIMP2)
     $ ,V4(NDIMP2),V5(NDIMP2),ENDELTA(NDIMP2),SCOLD(NFDIM,NDDIM)
     $ ,VOLD(NDIMP2),DB(NDIMP2,NDIMP2),WCHARM(NDDIM,NFDIM)
     $ ,DM(NDIMP2,NDIMP2),CRATE(NDIM,NDIM),RRATE(NDIM,NDIM)
     $ ,ETA(NDDIM),OPA(NDDIM),THOMSON(NDDIM),TAUTHOM(NDDIM)
     $ ,TNEW(NDDIM),XJCAPP(NFDIM),OPAC(NFDIM),SCNEW(NFDIM),DOPA(NFDIM)
     $ ,DETA(NFDIM),ETAC(NFDIM),EXPFAC(NFDIM),SIGMAKI(NFDIM,NDIM)
     $ ,PHI(NFLDIM),PWEIGHT(NFLDIM),DEPART(NDDIM,NDIM)
!     $ ,KODAT(MAXATOM),NFIRST(MAXATOM),NLAST(MAXATOM)
!     $ ,NCHARG(NDIM),MAINQN(NDIM),NOM(NDIM),ITNE(NDDIM),NFEDGE(NDIM)
!     $ ,MAINQN(NDIM),NOM(NDIM)
     $ ,ITNE(NDDIM),NFEDGE(NDIM)
     $ ,IWARN(NDDIM),LEVELPL(NDIM),MODHIST(MAXHIST)

      parameter (IPDIM=25,NBDIM=99)

!      character*8 :: agaunt

      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE
!      COMMON /COMIND/  INDNUP(MAXIND),INDLOW(MAXIND)
      COMMON /COMIND/
     $ XRED(MAXIND),XBLUE(MAXIND),XJLAPP(MAXIND)
     $ ,DETAL(MAXIND),OPAL(MAXIND),DOPAL(MAXIND)
     $ ,SCOLIND(MAXIND)
      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2
      integer   NGAMR(10),NGAML(10)
      real*8    AGAMR(10),AGAML(10)

      integer   XLBKG1,XLBKG2

      logical   LINE(MAXIND),NOTEMP,TPLOT,NODM,LBKG
      character MODHEAD*104,MODHOLD*104, CARD*80, LCARD*120
!      character LEVEL(NDIM)*10
      character*10 MAINPRO(NDDIM),MAINLEV(NDDIM)
!      character*10 ELEMENT(MAXATOM)
!      character*4 KEYCOL(NDIM,NDIM)
!      character*2 SYMBOL(MAXATOM)
      character*7 JOB
      integer,external :: time
      integer :: tstart

!      dimension WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)
      dimension xneclc(nddim)

      logical PROLIB

      REAL*8 :: CORMAX

      WRITE(*, *), 'STEAL: Entering STEAL, JOB = '//JOB

      IF(LBKG) PRINT*, 'STEAL: LINE BLANKETING SET TO TRUE'

      tstart = time()

C***  READING THE ATOMIC DATA FROM FILE DATOM

	CALL DATOM_M(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $               EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $               INDNUP,INDLOW,LASTIND,NATOM,
     $               ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $               NLAST,WAVARR,SIGARR,NFDIM)

C***  DECODING INPUT DATA ******************************************
      CALL DECSTE(LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,
     $     IFRRA,ITORA,IPRICC,IPRILC,LSEXPO,
     $     IFLUX,IDAT,LEVELPL,NDIM,IPLOTF,NEWWRC,
     $     NGAMR,NGAML,AGAMR,AGAML,DELTAC,LINE,MAXIND,NOTEMP,TPLOT,
     $     Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB,LBLANK)

C***  READING OF THE MODEL FILE ----------------------------------------
      IFL = 3; open(IFL, file='MODFILE', STATUS='OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,
     $             NDDIM,NPDIM,NFDIM,LBLANK)

      close(IFL)

      IFL = 3; open(IFL, file='POPNUM', STATUS='OLD')

c***  pop1 is dummy read because it will be overwritten below
      call readpop(ifl,T,popnum,pop2,pop3,pop1,rne,n,nd,modhead,jobnum)

      close(IFL)

c***  read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)	

      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             NCHARG,EDDARR,EDDI,NOM,WCHARM,
     $		   N,lastind,EINST,MODHEAD,JOBNUM)

c***  advance job-number counter
      JOBNUM = JOBNUM + 1

      CALL change (popnum,pop1,nd*n)

      CALL REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)
      if (lblank.lt.0) then
c***     the new blanketing table needs to be written to the model file
         IFL=3
         open (IFL,file='MODFILE',STATUS='UNKNOWN')
         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)
         CLOSE (ifl)
      endif

      DO 1 MG=1,10
      IF (NGAMR(MG).LE.JOBNUM) GAMMAR=AGAMR(MG)
      IF (NGAML(MG).LE.JOBNUM) GAMMAL=AGAML(MG)
    1 CONTINUE

C***  IN CASE OF JOBNUM=1 (STARTJOB), GAMMA'S ARE SET ZERO
      IF (JOBNUM .LE. 1) THEN

         GAMMAR = 0.0D0
         GAMMAL = 0.0D0
         NOTEMP = .TRUE.

      ENDIF
 
      IF (GAMMAR .EQ. 0. .AND. GAMMAL .EQ. 0. .AND. NOTEMP) THEN

C***  CALCULATION OF POPNUMBERS WITH THE USUAL (LINEAR) RATE EQUATION.
C***  THIS IS ECONOMIC AND PROVIDES THE ORIGINAL RATES FOR A POSSIBLE
C***  PRINTOUT BY SUBR. PRIRAT
      
      CALL POPZERO(T,RNE,POPNUM,DEPART,ENTOT,ITNE,N,ENLTE,
     $             WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     $             XLAMBDA,FWEIGHT,XJCARR,NF,XJL,IFRRA,ITORA,ALPHA,
     $             SEXPO,AGAUNT,MODHEAD,MODHOLD,JOBNUM,
     $             LASTIND,ND,LSRAT,CRATE,RRATE,RATCO,
     $             SIGMAKI,ALTESUM,COCO,KEYCOL,NOM,NATOM,
     $             KODAT,NFIRST,NLAST,WAVARR,SIGARR)

C***  PROVIDING THE TEMPERATURE STRATIFICATION T(R) FOR A POSSIBLE OUTPUT BY
C***  SUBR.S PRITAU, PLOTT

      TNEW(1 : ND) = T(1 : ND)

      ELSE

C***  RATE EQUATION WITH APPROXIMATE RADIATION TRANSFER (SCHARMER-METHOD)
C***  IN THIS BRANCH, PRIRAT MAY ONLY SHOW THE NETTO RATES
C***  CALCULATION OF NEW POPULATION NUMBERS, EL. DENSITY AND DEPARTURE COEFF.

      CALL LINPOP(T,RNE,ENTOT,ITNE,POPNUM,DEPART,POP1,
     $   N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EN,ENDELTA,EINST,LEVEL,
     $   XLAMBDA,FWEIGHT,XJCARR,NF,XJL,WCHARM,XJCAPP,SCOLD,XJLAPP,
     $   DM,DB,V1,V2,V4,V5,VOLD,GAMMAL,EPSILON,
     $   TNEW,NOTEMP,NODM,IADR19,MAXADR,
     $   DELTAC,GAMMAR,IPRICC,IPRILC,MODHEAD,JOBNUM,IFRRA,ITORA,
     $   RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $   VELO,GRADI,VDOP,PHI,PWEIGHT,SCOLIND,INDNUP,INDLOW,LASTIND,
     $   OPAC,SCNEW,DOPA,DETA,OPAL,DOPAL,DETAL,SIGMAKI,
     $   ND,LSRAT,CRATE,RRATE,RATCO,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,
     $   LINE,ALTESUM,ETAC,NFEDGE,EXPFAC,NOM,NATOM,KODAT,NFIRST,
     $   NLAST,WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,JOBMAX)

      ENDIF
 
C***  REDUCED CORRECTIONS, IF OPTION IS SET
      IF (REDUCE .NE. 1.0D0 .AND. JOBNUM .GT. 1) CALL REDCOR (POPNUM,POP1,ND,N,RNE,NCHARG,REDUCE)
 
      IF (LSPOP.GT.0) CALL PRIPOP (LSPOP,WEIGHT,NCHARG,NOM,TNEW,NOTEMP,
     $             ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)
 
      IF (JOBNUM .GT. 1) CALL PRICORR(POPNUM,POP1,LEVEL,N,ND,MODHEAD,LSPOP,CORMAX,
     $                                NCHARG,RNE,JOBNUM,REDUCE,GAMMAL,GAMMAR,
!     $                                T,TNEW,NOTEMP,DELTAC,ELEMENT,NATOM,NFIRST,NLAST,LBKG)
     $                                T,TNEW,NOTEMP,DELTAC,ELEMENT,NATOM,NFIRST,NLAST)

!      IF (JOBNUM .GT. 1) LAMBDA_ITER = LAMBDA_ITER + 1
 
      IF (LSEXPO .GT. 0.AND.JOBNUM.GT.3) CALL PRIEXPO(POPNUM,POP1,POP2,LEVEL,N,ND,MODHEAD,JOBNUM,LSEXPO)
 
C***  UPDATING THE MODEL FILE
      ifl = 2; OPEN(ifl, FILE = 'POPNUM', STATUS = 'unknown')

      XNECLC(1 : ND) = RNE(1 : ND) * ENTOT(1 : ND)

      call writpop(ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum); CLOSE(ifl)

c***  if a new LB table is read then store a new model file
      if (lblank.lt.0) then
      IFL=3
      open (IFL,file='MODFILE',STATUS='UNKNOWN')
      CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
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
      CALL STHIST(LCARD,GAMMAL,GAMMAR,NOTEMP,DELTAC,MODHEAD,JOBNUM,CORMAX,REDUCE,MODHOLD,time()-tstart)
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

            CALL PRIPOP(LSPOP,WEIGHT,NCHARG,NOM,TNEW,NOTEMP,
     $                  ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD)

            print*, 'steal: ', NF; stop

            CALL PRITAU(MODHEAD,JOBNUM,RSTAR,ND,RADIUS,RNE,ENTOT,TNEW,
     $                  POPNUM,N,EN,LEVEL,NCHARG,WEIGHT,ELEVEL,
     $                  EION,EINST,ALPHA,SEXPO,AGAUNT,NOM,XLAMBDA,
     $                  FWEIGHT,TAUTHOM,TAUROSS,WAVARR,SIGARR,NF)

            CALL PRIH(LPRIH,ND,RADIUS,HTOT,TEFFE,
     $                TNEW,TAUROSS,JOBNUM,MODHEAD)
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
     $            KODAT,AKEY,NF,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
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
      IF (TPLOT) CALL PLOTT (ND,RADIUS,TNEW,MODHEAD,JOBNUM)
 
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

      END MODULE
