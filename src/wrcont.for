      MODULE MOD_WRCONT

      CONTAINS

      SUBROUTINE WRCONT(job)

      use MOD_AMBIPOLAR
      use MOD_CALCH
      use MOD_COOP
      use MOD_DATOM
      use MOD_DECON
      use MOD_DECSTAR
      use MOD_DIFFUS
      use MOD_ELIMIN
      use MOD_extrxjc
      use MOD_FORMATS
      use MOD_PHNU
      use MOD_PRIBLA
      use MOD_PRIFLUX
      use MOD_PRIGH
      use MOD_PRIINT
      use MOD_PRIOPA
      use MOD_PRIV
      use MOD_READMOD
      use MOD_READPOP
      use MOD_READRAD
      use MOD_REBLANK
      use MOD_STORXJC
      use MOD_TICTOC
      use MOD_WRITMOD
      use MOD_WRITRADC
      use ABUNDANCES
      USE CONSTANTS

      use common_block
      use file_operations
      use vardatom
      use varhminus
      use varsteal

      use local_operator
 
!***  THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
!***  CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
!***  IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
!***  AND THE FREQUENCY GRID (FILE FGRID)
!***  PRESENT VERSION: MODEL ATMOSPHERE OF HELIUM (CODE NR. "1") WITH
!***                                       HYDROGEN         "2"
!***  FOR IMPLEMENTATION OF ADDITIONAL ELEMENTS:
!***  MODIFY SUBROUTINES  "DATOM", "DECSTAR"
!***  INSERT CORRESPONDING ATOMIC DATA INTO SUBR. "COLLI", "PHOTOCS"
      IMPLICIT REAL*8(A - H, O - Z)
 
      real*8, ALLOCATABLE :: DUMMY2(:, :)

!***  CHANGES BY MARGIT HABERREITER, 20 MAY, 2002
!***  LEVLOW NEEDS TO BE DEFINED, AS IT IS USED AS A KEYWORD TO SELECT THE 
!***  ELEMENT AND LEVEL TO READ THE CONTINUUM OPACITIES FROM AN INPUT TABLE
CMH   LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH   XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF

!     THIS PROGRAM IS TO SOLVE THE CONTINUOUS RADIATION TRANSFER
!     WITH GIVEN POPULATION NUMBERS

      COMMON // U(NDDIM,NPDIM),VL(NPDIM),HNU(NDDIM),VJL(NPDIM,NDDIM)

      parameter (IPDIM = 25, NBDIM = 99)
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)

      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2 
      LOGICAL LBKG
      INTEGER XLBKG1,XLBKG2 
      CHARACTER MODHEAD*104,CARD*80, LCARD*100
      CHARACTER*10 MAINPRO(NDDIM),MAINLEV(NDDIM)

      CHARACTER*7 JOB

      integer :: timer

      real*8 :: amu

      LOGICAL LDUMMY1, LDUMMY2

      REAL*8, ALLOCATABLE, DIMENSION(:, :) :: EDDI_OLD

      REAL*8 ::  DEDDI1, DEDDI2, DEDDI3

      INTEGER :: DEDDI1_LOC, DEDDI2_LOC, DEDDI3_LOC

      DATA AMU /1.660531d-24/

      print*, 'Entering wrcont, JOB = '//JOB

      CALL DATOM(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $           EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $           INDNUP,INDLOW,LASTIND,NATOM,
     $           ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $           NLAST,WAVARR,SIGARR,NFDIM)

      CALL DECSTAR(MODHEAD,FM,RSTAR,teff,glog,xmass,VDOP,LDUMMY1,LBKG,XLBKG1,XLBKG2,
     $             LDUMMY2,NATOM,ABXYZ,KODAT,IDAT,LBLANK,ATMEAN,amu)

!***  DECODING INPUT OPTIONS *******************************************
      CALL DECON(LSOPA,LSINT,IFLUX,JOBMAX,LPRIH,LPHNU,LPRIV,TEFF,LBLANK)

!***  READING OF THE MODEL FILE ----------------------------------------
      IFL = 3

      open(IFL, file = 'MODFILE', STATUS = 'OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      if (allocated(dummy2)) deallocate(dummy2); allocate(dummy2(NF, N))

      IF (JOBNUM .GE. 1000) JOBNUM = JOBNUM - 100

      IFL = 3

      open(IFL, file = 'POPNUM', STATUS = 'OLD')

      call readpop(ifl, T, popnum, pop1, pop2, pop3, rne, n, nd, modhead, jobnum)

      close(ifl)

!     read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage) 
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             NCHARG,EDDARR,EDDI,NOM,WCHARM,
     $             N,lastind,EINST,MODHEAD,JOBNUM)

      JOBNUM = JOBNUM + 1

!     initialize frequency integrated quantities
      TOTIN =        0.0d0
      TOTOUT =       0.0d0

      GTOT(1 : ND) = 0.0d0
      ETOT(1 : ND) = 0.0d0
      XTOT(1 : ND) = 0.0d0
      HTOT(1 : ND) = 0.0d0

!***  the blanketing table is read by routine READMOD
!***  if lblank.gt.0 then read a new table from the file LIBLANK
      CALL REBLANK(LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

      if (lblank.lt.0) then
!***     the new blanketing table needs to be written to the model file
         IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'UNKNOWN')

         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,
     $                XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $                ABXYZ,NATOM,MODHEAD,JOBNUM) 

         CLOSE(ifl)

      endif

      IF (abs(LBLANK).EQ.2) CALL PRIBLA(LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,SCAFAC,ABSFAC)

      call TIC(timer)

!***  SOLUTION OF THE TRANSFER EQUATION FOR EACH FREQUENCY-POINT ********

      ALLOCATE(EDDI_OLD(3, ND)); EDDI_OLD = EDDI(1 : 3, 1 : ND)

      FRQS: DO K = 1, NF

!       now extract XJC and EDDI for the frequency K
        CALL EXTRXJC(XJCARR, XJC, EDDARR, EDDI, nd, nf, K)

        CALL COOP(XLAMBDA(K),ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $            OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $            N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $            ALPHA,SEXPO,AGAUNT,0,DUMMY2,
     $            WAVARR(1 : N, 1 : NF), SIGARR(1 : N, 1 : NF),
     $            LBKG,XLBKG1,XLBKG2,NF)

!        print*, 'wrcont: K = ', K, ' out of ', NF

        CALL DIFFUS(XLAMBDA(K),T,RADIUS,ND,BCORE,DBDR)

        if (LSOPA .GT. 0) CALL PRIOPA(XLAMBDA(K),K,ND,LSOPA,RADIUS,
     $                                OPA,ETA,THOMSON,IWARN,MAINPRO,
     $                                MAINLEV,JOBNUM,MODHEAD)

        if (any(isnan(eddi(1, :)))) then

            do i = 1, ND

               write(*, '(i4,2x,e15.7,2x,i4,3(2x,e15.7))') 
     $         k, fweight(k), i, etot(i), eddi(1, i), xjc(i)

            enddo

            stop 'eddi nan before elimin'

        endif

!        print*, 'before elimin'

        CALL ELIMIN(XLAMBDA(K),EMFLUX(K),FLUXIN,U,Z,XJC,RADIUS,
     $              P,BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

!        print*, 'after elimin'

        !***  INTEGRATION OF THE TOTAL INCIDENT AND EMERGENT FLUX
        TOTIN=TOTIN+FLUXIN*FWEIGHT(K)
        TOTOUT=TOTOUT+EMFLUX(K)*FWEIGHT(K)

        CALL CALCH(ND,NP,NPDIM,OPA,Z,P,U,VL,VJL,RADIUS,HNU,FLUXIN,EMFLUX(K))

        XTOT(1 : ND) = XTOT(1 : ND) + XJC(1 : ND) * FWEIGHT(K)

        if (any(isnan(eddi(1, :)))) then

            do i = 1, ND

               write(*, '(i4,2x,e15.7,2x,i4,3(2x,e15.7))') 
     $         k, fweight(k), i, etot(i), eddi(1, i), xjc(i)

            enddo

            stop 'eddi nan before'

        endif

        ETOT(1 : ND) = ETOT(1 : ND) + EDDI(1, 1 : ND) * XJC(1 : ND) * FWEIGHT(K)

        if (any(isnan(etot))) then

            do i = 1, ND

               write(*, '(i4,2x,e15.7,2x,i4,3(2x,e15.7))') 
     $         k, fweight(k), i, etot(i), eddi(1, i), xjc(i)

            enddo

            stop 'nan after'

        endif

        GTOT(1 : ND) = GTOT(1 : ND) + OPA(1 : ND) * HNU(1 : ND) * FWEIGHT(K)

        HTOT(1 : ND) = HTOT(1 : ND) + HNU(1 : ND) * FWEIGHT(K)

!       XJC and EDDI are stored for later write to file RADIOC
        call storxjc(XJCARR, XJC, EDDARR, EDDI, nd, nf, K)

!       PRINTOUT OF FREQUENCY DEPENDEND VARIABLES
        IF (LPRIV.GT.0.AND.(K.LT.68 .OR. K.EQ.111)) CALL PRIV(K,XLAMBDA(K),LPRIV,ND,NP,RADIUS,Z,VJL,RSTAR)

        IF (LPHNU.GT.0) CALL PHNU(K,XLAMBDA(K),LPHNU,ND,RADIUS,HNU,RSTAR)

!       ABORT, IF NOT SUFFICIENT TIME LEFT FOR NEXT FREQUENCY POINT
        LASTK = K

      ENDDO FRQS

      print*, 'maxmin of dEDDI(1, :) = ', 
     & maxval(abs(EDDI(1, 1 : ND) / EDDI_OLD(1, :))),
     & minval(abs(EDDI(1, 1 : ND) / EDDI_OLD(1, :)))

      print*, 'maxmin of dEDDI(2, :) = ', 
     & maxval(abs(EDDI(2, 1 : ND) / EDDI_OLD(2, :))),
     & minval(abs(EDDI(2, 1 : ND) / EDDI_OLD(2, :)))

      print*, 'maxmin of dEDDI(3, :) = ', 
     & maxval(abs(EDDI(3, 1 : ND) / EDDI_OLD(3, :))),
     & minval(abs(EDDI(3, 1 : ND) / EDDI_OLD(3, :)))

!=============================================================================
!Rinat Tagirov

      DEDDI1 =     MAXVAL(ABS(EDDI(1, 1 : ND) / EDDI_OLD(1, :)) - 1.0D0)

      DEDDI1_LOC = MAXLOC(ABS(EDDI(1, 1 : ND) / EDDI_OLD(1, :)) - 1.0D0, 1)

      DEDDI2 =     MAXVAL(ABS(EDDI(2, 1 : ND) / EDDI_OLD(2, :)) - 1.0D0)

      DEDDI2_LOC = MAXLOC(ABS(EDDI(2, 1 : ND) / EDDI_OLD(2, :)) - 1.0D0, 1)

      DEDDI3 =     MAXVAL(ABS(EDDI(3, 1 : ND) / EDDI_OLD(3, :)) - 1.0D0)

      DEDDI3_LOC = MAXLOC(ABS(EDDI(3, 1 : ND) / EDDI_OLD(3, :)) - 1.0D0, 1)

      IF (LAMBDA_ITER .NE. 0) THEN

         CALL OPEN_TO_APPEND(1836, EDDI_FILE)

         WRITE(1836, '(I5,3(2x,E15.7,2x,I3))') LAMBDA_ITER, DEDDI1, DEDDI1_LOC, DEDDI2, DEDDI2_LOC, DEDDI3, DEDDI3_LOC

         CLOSE(1836)

      ENDIF

!=============================================================================================================================

      call AMBIPOLAR(ND, N, T, ENTOT, RNE, LEVEL, RADIUS, POPNUM, RSTAR, timer)

!     store the new continuum radiation field
      IF (LASTK.EQ.NF) ETOT(1:ND)=ETOT(1:ND)/XTOT(1:ND)

      call writradc(xjcarr,xjc,eddarr,eddi,emflux,totin,totout,HTOT,GTOT,XTOT,ETOT,wcharm,nd,nf,MODHEAD,JOBNUM)

!     if a new LB table is read then store a new model file
      if (lblank.lt.0) then

         IFL = 3

         open(IFL, file = 'MODFILE', STATUS='UNKNOWN')

         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,
     $                XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)

         close(ifl)

      endif

      open (7,file='MODHIST',status='old')

      do
        read (7,'(A80)',end=11) card
      enddo
 11   continue

!***  UPDATING THE MODEL HISTORY
      write(LCARD,88) JOBNUM,'. WRCONT   LASTK=',LASTK,TOC(TIMER),' sec'
   88 FORMAT (1H/,I3,A,I5,' COMPLETE  -  run time all K: ',i10,A)
      write (7,'(A100)') LCARD
      close (7)

      !***  PRINTOUTS
      IF (LSINT.GT.0 .AND. LASTK .EQ. NF)
     $   CALL PRIINT (XJC,XJCARR,EDDI,EDDARR,RADIUS,ND,XLAMBDA,NF,
     $                LSINT,JOBNUM,MODHEAD)
      IF (IFLUX.GT.0 .AND. LASTK .EQ. NF)
     $   CALL PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,JOBNUM,
     $                 FWEIGHT,MODHEAD,AKEY )
      ! IF (LPRIH.GT.0 .AND. LASTK .EQ. NF)
      IF (LPRIH.GT.0)
     $   CALL PRIGH (LPRIH,ND,RADIUS,HTOT,GTOT,ETOT,TEFF,ENTOT,RNE,
     $               RSTAR,T,VELO,GRADI,ATMASS,ABXYZ,NATOM)






      !***  ROUTING OF SUBSEQUENT JOBS
      IF (LASTK .EQ. NF) THEN
         WRITE (6,*) ' ALL FREQUENCY POINTS COMPLETED'
         IF (JOBNUM .GE. JOBMAX) THEN
            JOB='exit'
            PRINT *,' MAX. NUMBER OF JOBS EXCEEDED, JOB='//JOB
            ELSE
            JOB='newline'
            PRINT *,' REPEAT JOB TO BE ROUTED, JOB='//JOB
            ENDIF
         ELSE
         JOB='wrcont'
         PRINT *,LASTK,' OF ',NF,' FREQUENCY POINTS COMPLETED'
         PRINT *,' NEW WRCONT-JOB TO BE ROUTED, JOB='//JOB
         ENDIF

      RETURN

      END subroutine

      end module
