      MODULE MOD_WRCONT

      CONTAINS

      SUBROUTINE WRCONT (job)

      use MOD_AMBIPOLAR
      use MOD_CALCH
      use MOD_COOP_M
      use MOD_CPPOPNUM
      use MOD_DATOM_M
      use MOD_DECON
      use MOD_DECSTAR_M
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

      USE MOD_CALCLAMBDAS
      USE CONSTANTS

      USE COMMON_BLOCK
      USE FILE_OPERATIONS

      use VARDATOM
 
C***********************************************************************
C***  THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
C***  CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
C***  IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
C***  AND THE FREQUENCY GRID (FILE FGRID)
C***  PRESENT VERSION: MODEL ATMOSPHERE OF HELIUM (CODE NR. "1") WITH
C***                                       HYDROGEN         "2"
C***  FOR IMPLEMENTATION OF ADDITIONAL ELEMENTS:
C***                          MODIFY SUBROUTINES  "DATOM", "DECSTAR"
C***  INSERT CORRESPONDING ATOMIC DATA INTO SUBR. "COLLI", "PHOTOCS"
C***********************************************************************
      IMPLICIT REAL*8(A - H, O - Z)
 
      real*8, ALLOCATABLE :: DUMMY2(:,:)

C******************************************************************************
C***  CHANGES BY MARGIT HABERREITER, 20 MAY, 2002
C***  LEVLOW NEEDS TO BE DEFINED, AS IT IS USED AS A KEYWORD TO SELECT THE 
C***  ELEMENT AND LEVEL TO READ THE CONTINUUM OPACITIES FROM AN INPUT TABLE
C      CHARACTER*10 LEVLOW
CMH   DIMENSION WAVARR(NDIM,97),SIGARR(NDIM,97),WLTH(NDIM,97)
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF

!      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)

c     INTEGER XLBKG1,XLBKG2
c     LOGICAL LBKG
C******************************************************************************
c      SUBROUTINE WRCONT (JOB)
C*******************************************************************************
***  THIS PROGRAM IS TO SOLVE THE CONTINUOUS RADIATION TRANSFER
C***  WITH GIVEN POPULATION NUMBERS
C*******************************************************************************
      COMMON // RADIUS(NDDIM),ENTOT(NDDIM),T(NDDIM)
     $ ,XJC(NDDIM),XJCARR(NDDIM,NFDIM),XJL(NDDIM,MAXIND)
     $ ,EDDI(3,NDDIM),EDDARR(3,NDDIM,NFDIM),TAUROSS(NDDIM)
     $ ,RNE(NDDIM),VELO(NDDIM),GRADI(NDDIM),AINCRIT(NDDIM)
     $ ,XLAMBDA(NFDIM),FWEIGHT(NFDIM),EMFLUX(NFDIM),AKEY(NFDIM)
!     $ ,WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
!     $ ,EINST(NDIM,NDIM),ALPHA(NDIM),SEXPO(NDIM)
     $ ,ENLTE(NDIM)
!     $ ,COCO(NDIM,NDIM,4),ALTESUM(4,NDIM)
     $ ,P(NPDIM),Z(NDDIM,NPDIM),POPNUM(NDDIM,NDIM)
!     $ ,ATMASS(MAXATOM),STAGE(MAXATOM),AGAUNT(NDIM)
     $ ,HTOT(NDDIM),GTOT(NDDIM),XTOT(NDDIM),ETOT(NDDIM)
     $ ,POP1(NDDIM,NDIM),POP2(NDDIM,NDIM),POP3(NDDIM,NDIM)
c***  up to here the common is identical to that of steal
c      COMMON // WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
c     $ ,MAINQN(NDIM),EINST(NDIM,NDIM),ALPHA(NDIM),SEXPO(NDIM)
c     $ ,IGAUNT(NDIM),ALTESUM(4,NDIM),NOM(NDIM)
c     $ ,COCO(NDIM,NDIM,4)
c     $ ,ABXYZ(MAXATOM),KODAT(MAXATOM),ATMASS(MAXATOM),STAGE(MAXATOM)
c     $ ,NFIRST(MAXATOM),NLAST(MAXATOM)
c     $ ,OPA(NDDIM),ETA(NDDIM),THOMSON(NDDIM),IWARN(NDDIM)
c     $ ,VELO(NDDIM),GRADI(NDDIM)
c     $ ,RADIUS(NDDIM),ENTOT(NDDIM),T(NDDIM),RNE(NDDIM),XJC(NDDIM)
c     $ ,EDDI(3,NDDIM),IADR7(NDDIM)
     $ ,ETA(NDDIM),OPA(NDDIM),THOMSON(NDDIM),TAUTHOM(NDDIM)
C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY !
     $ ,A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM)
     $ ,BX(NPDIM,NPDIM,NDDIM),WX(NPDIM,NDDIM)
c     $ ,P(NPDIM),A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM)
c     $ ,BX(NPDIM,NPDIM,NDDIM),WX(NPDIM,NDDIM)
c     $ ,XLAMBDA(NFDIM),EMFLUX(NFDIM),FWEIGHT(NFDIM),KEY(NFDIM)
c     $ ,U(NDDIM,NPDIM),Z(NDDIM,NPDIM),VL(NPDIM),HNU(NDDIM),HTOT(NDDIM)
c     $ ,POPNUM(NDDIM,NDIM),VJL(NPDIM,NDDIM),GTOT(NDDIM),XTOT(NDDIM)
     $ ,U(NDDIM,NPDIM),VL(NPDIM),HNU(NDDIM),VJL(NPDIM,NDDIM)
c*** the integers in the main common are identical to that in steal
!     $ ,KODAT(MAXATOM),NFIRST(MAXATOM),NLAST(MAXATOM)
!     $ ,NCHARG(NDIM),MAINQN(NDIM),NOM(NDIM)
     $ ,ITNE(NDDIM),NFEDGE(NDIM)
     $ ,IWARN(NDDIM),LEVELPL(NDIM),MODHIST(MAXHIST)
c***  in steal the COMIND common is larger
!      COMMON /COMIND/  INDNUP(MAXIND),INDLOW(MAXIND)

c     $ ,INDNUP(MAXIND),INDLOW(MAXIND),      ETOT(NDDIM)
c     $ ,MODHIST(MAXHIST),IADR3(MAXADR),NCHARG(NDIM)
      parameter (IPDIM = 25, NBDIM = 99)
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)
C23456789 12345678901234567890123456789012345678901234567890123456789012
      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2 
      LOGICAL LBKG
      INTEGER XLBKG1,XLBKG2 
      CHARACTER MODHEAD*104,CARD*80, LCARD*100
      CHARACTER*10 MAINPRO(NDDIM),MAINLEV(NDDIM)
!      CHARACTER LEVEL(NDIM)*10
!      CHARACTER*10 ELEMENT(MAXATOM)
!      CHARACTER*4 KEYCOL(NDIM,NDIM)
!      CHARACTER*2 SYMBOL(MAXATOM)

      CHARACTER*7 JOB

!      character*8 :: agaunt

      LOGICAL NOTEMP

      integer :: timer

      REAL*8, DIMENSION(NDDIM, NFDIM) :: WCHARM

      LOGICAL LDUMMY1, LDUMMY2, LDUMMY3

      REAL*8, ALLOCATABLE, DIMENSION(:, :) :: EDDI_OLD

!      REAL*8, ALLOCATABLE, DIMENSION(:) :: OpticalDepth

!      REAL*8, ALLOCATABLE, DIMENSION(:) :: LambdaOperatorDiag

      REAL*8 ::  DEDDI1, DEDDI2, DEDDI3

      INTEGER :: DEDDI1_LOC, DEDDI2_LOC, DEDDI3_LOC

C***  write output to a file
      print *, 'wrcont : Entering WRCONT  , JOB='//JOB
!      open (6,file='wrcont.out',status='unknown')
C**** Changes by Margit Haberreiter
CMH   CALL       DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
       CALL      DATOM_M(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                   EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $                   INDNUP,INDLOW,LASTIND,NATOM,
     $                   ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $                   NLAST,WAVARR,SIGARR,NFDIM)

C***  MASS STORAGE ON FILE 7 FOR THE FEAUTRIER MATRICES
c      IERR=1
C      CALL OPENMS(7,IADR7,NDDIM,0,IERR)
c      IF (IERR.LT.0) THEN
c         PRINT *, ' IERR = ',IERR
c         STOP 'OPMS7'
c         ENDIF
      CALL DECSTAR_M(MODHEAD,FM,RSTAR,VDOP,RMAX,LDUMMY1,LDUMMY2,LBKG,XLBKG1,XLBKG2,
     $               LDUMMY3,NATOM,ABXYZ,KODAT,IDAT,LBLANK,ATMEAN)
C***  DECODING INPUT OPTIONS *******************************************
      CALL       DECON(LSOPA,LSINT,IFLUX,JOBMAX,
     $                 LPRIH,LPHNU,LPRIV,TEFF,NOTEMP,LBLANK)

C***  READING OF THE MODEL FILE ****************************************
c      CALL       RMODCON (ND,NDDIM,RADIUS,NP,NPDIM,P,Z,ENTOT,T,RNE,NF,
c     $             NFDIM,MODHIST,MAXHIST,LAST,ALTESUM,XLAMBDA,
c     $             IADR3,MAXADR,
c     $             FWEIGHT,AKEY,POPNUM,RSTAR,MODHEAD,JOBNUM,NEXTK,N,
c     $             VELO,GRADI,ABXYZ,NATOM)

C***  READING OF THE MODEL FILE ----------------------------------------
      IFL=3
      open (IFL,file='MODFILE',STATUS='OLD')
      CALL READMOD       (IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                    GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                    ABXYZ,NATOM,MODHEAD,JOBNUM,
     $                    NDDIM,NPDIM,NFDIM,NEXTK,LBLANK)
      close (IFL)

      print*, 'wrcont: ', NF, NFDIM

      stop

      IF (JOBNUM .GE. 1000) JOBNUM = JOBNUM - 100

      IFL=3
      open (IFL,file='POPNUM',STATUS='OLD')
      call readpop (ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,
     $                    jobnum)
      close (ifl)
c***  read the radiation field from files RADIOC and RADIOL
c          (pop1 is used as dummy storage) 
      CALL READRAD    (NF,ND,POP1 ,XJCARR,XJC,XJL,
     $                 HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $                 NCHARG,EDDARR,EDDI,NOM,WCHARM,
     $                 N,MAXIND,EINST,NDIM,MODHEAD,JOBNUM)

      JOBNUM = JOBNUM + 1

c      missing:     ALTESUM,
C***   initialize frequency integrated quantities
      TOTIN=.0
      TOTOUT=.0
      GTOT(1:ND)=0.
      ETOT(1:ND)=0.
      XTOT(1:ND)=0.
      HTOT(1:ND)=0.
      IF (NEXTK .ne. 1) THEN
c        there was once the option to continue the frequency loop...
         PRINT *,' NEXTK= ',NEXTK
         write (6,*) 'this operation is no longer active'
         pause
         stop
      ENDIF
 
c***  the blanketing table is read by routine READMOD
c***  if lblank.gt.0 then read a new table from the file LIBLANK
      CALL REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

      if (lblank.lt.0) then
c***     the new blanketing table needs to be written to the model file
         IFL=3
         open (IFL,file='MODFILE',STATUS='UNKNOWN')
         CALL WRITMOD       (IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                    GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                    ABXYZ,NATOM,MODHEAD,JOBNUM) 
         CLOSE (ifl)
      endif
      IF (abs(LBLANK).EQ.2) 
     $  CALL PRIBLA (LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,
     $                   SCAFAC,ABSFAC)
      call TIC(timer)

C***  SOLUTION OF THE TRANSFER EQUATION FOR EACH FREQUENCY-POINT ********

!      ALLOCATE(LambdaOperatorDiag(ND))
!      ALLOCATE(OpticalDepth(ND))

!      OPEN (UNIT = 10, FILE = 'continuum_lambda_operator_diag.out', 
!     $ FORM = 'FORMATTED')

!      WRITE(10, 99999)

!99999 FORMAT('Wavelength',10x,'ID',10x,'Height',10x,
!     $ 'Opacity',10x,'Emissivity',10x,
!     $ 'Optical Depth',10x,'LambdaOperDiag')
 
      ALLOCATE(EDDI_OLD(3, ND));
      EDDI_OLD = EDDI(1 : 3, 1 : ND)

      FRQS: DO K = NEXTK, NF

        !*** now extract XJC and EDDI for the frequency K

        CALL EXTRXJC (XJCARR,XJC,EDDARR,EDDI,nd,nf,K)

        CALL COOP_M(XLAMBDA(K),ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $              OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $              N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $              ALPHA,SEXPO,AGAUNT,0,DUMMY2,WAVARR,SIGARR,
     $              LBKG,XLBKG1,XLBKG2,NF)

        CALL DIFFUS (XLAMBDA(K),T,RADIUS,ND,BCORE,DBDR)

        if (LSOPA.GT.0) then

          CALL PRIOPA (XLAMBDA(K),K,ND,LSOPA,RADIUS,
     $      OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)

        endif

!        OpticalDepth(1) = 0.0D0

!        DO L = 2, ND

!           OpticalDepth(L) = OpticalDepth(L - 1) + 
!     $ ((OPA(L) + OPA(L - 1)) / 2.) * ((Height(L - 1) - Height(L)) / 
!     $ SolarRadiusKM)

!        ENDDO

!        OpticalDepth(1) = OpticalDepth(2)

!        LambdaOperatorDiag(:) = CALCLAMBDAS(OPA, RADIUS, EDDI, ND) ! Diagonal of continuum Lamda-operator

!        DO L = 1, ND

!          WRITE (10,'(E15.7,8x,I2,8x,E15.7,8x,E15.7,8x,E15.7,8x,E15.7,
!     $  8x,E15.7)')
!     $  XLAMBDA(K), L, Height(L), OPA(L), ETA(L),
!     $  OpticalDepth(L), LambdaOperatorDiag(L)

!        ENDDO

        CALL ELIMIN (XLAMBDA(K),EMFLUX(K),FLUXIN,U,Z,
     $          A,B,C,W,BX,WX,XJC,RADIUS,P,BCORE,DBDR,
     $                       OPA,ETA,THOMSON,EDDI,ND,NP,NPDIM)
        !***  INTEGRATION OF THE TOTAL INCIDENT AND EMERGENT FLUX
        TOTIN=TOTIN+FLUXIN*FWEIGHT(K)
        TOTOUT=TOTOUT+EMFLUX(K)*FWEIGHT(K)

        CALL CALCH(ND,NP,NPDIM,OPA,Z,P,U,VL,VJL,RADIUS,HNU,FLUXIN,
     $ EMFLUX(K))

        XTOT(:ND) = XTOT(:ND)+XJC(:ND)*FWEIGHT(K)
        ETOT(:ND) = ETOT(:ND)+EDDI(1,:ND)*XJC(:ND)*FWEIGHT(K)
        GTOT(:ND) = GTOT(:ND)+OPA(:ND)*HNU(:ND)*FWEIGHT(K)
        HTOT(:ND) = HTOT(:ND)+HNU(:ND)*FWEIGHT(K)

!      PRINT*, K, XLAMBDA(K), 'FREQUENCY CYCLE IN WRCONT.FOR IS GOING ON'

        !***  XJC and EDDI are stored for later write to file RADIOC
        call storxjc (XJCARR,XJC,EDDARR,EDDI,nd,nf,K)
        !***  WRITING ON THE MODEL FILE
        ! WRITE (CNAME,FMT_KEY) 'XJC ',K
        !! 3 FORMAT(A4,I4)
        ! CALL WRITMS(3,XJC,ND,NAME,-1,IDUMMY,IERR)
        ! WRITE (CNAME,3) 'EDDI',K
        ! CALL WRITMS (3,EDDI,3*ND,NAME,-1,IDUMMY,IERR)

        !***  PRINTOUT OF FREQUENCY DEPENDEND VARIABLES
        IF (LPRIV.GT.0.AND.(K.LT.68 .OR. K.EQ.111))
     $    CALL PRIV (K,XLAMBDA(K),LPRIV,ND,NP,RADIUS,Z,VJL,RSTAR)
        IF (LPHNU.GT.0)
     $    CALL PHNU (K,XLAMBDA(K),LPHNU,ND,RADIUS,HNU,RSTAR)

        !***  ABORT, IF NOT SUFFICIENT TIME LEFT FOR NEXT FREQUENCY POINT
        LASTK=K

      ENDDO FRQS
 
!      DEALLOCATE(LambdaOperatorDiag)
!      DEALLOCATE(OpticalDepth)

!      CLOSE(10)
!      CLOSE(110)

      print *,'maxmin of dEDDI(1,:) = ', 
     &  maxval(abs(EDDI(1,1:ND)/EDDI_OLD(1,:))),
     &  minval(abs(EDDI(1,1:ND)/EDDI_OLD(1,:)))
      print *,'maxmin of dEDDI(2,:) = ', 
     &  maxval(abs(EDDI(2,1:ND)/EDDI_OLD(2,:))),
     &  minval(abs(EDDI(2,1:ND)/EDDI_OLD(2,:)))
      print *,'maxmin of dEDDI(3,:) = ', 
     &  maxval(abs(EDDI(3,1:ND)/EDDI_OLD(3,:))),
     &  minval(abs(EDDI(3,1:ND)/EDDI_OLD(3,:)))

!============================================================================================================================
!RT

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

      call AMBIPOLAR(ND,N,T,ENTOT,RNE,LEVEL,RADIUS,POPNUM, RSTAR,timer)
      !***  UPDATING THE MODEL FILE
      ! CALL WRITMS (3,MODHIST,MAXHIST,7HMODHIST,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,EMFLUX,NF,6HEMFLUX,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,TOTIN,1,5HTOTIN,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,TOTOUT,1,6HTOTOUT,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,HTOT,ND,4HHTOT,-1,IDUMMY,IERR)
      !***  store the new continuum radiation field
      IF (LASTK.EQ.NF) ETOT(1:ND)=ETOT(1:ND)/XTOT(1:ND)
      call writradc (xjcarr,xjc,eddarr,eddi,emflux,totin,totout,
     $               HTOT,GTOT,XTOT,ETOT,wcharm,nd,nf,MODHEAD,JOBNUM)
      !***  if a new LB table is read then store a new model file
      if (lblank.lt.0) then
      IFL=3
      open (IFL,file='MODFILE',STATUS='UNKNOWN')
      CALL WRITMOD       (IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                    GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                    ABXYZ,NATOM,MODHEAD,JOBNUM)
      CLOSE (ifl)
      endif

C***  UPDATING THE MODEL HISTORY
c      ENCODE (24,8,MODHIST(LAST+1)) JOBNUM,LASTK
c      MODHIST(2)=JOBNUM
c     MODHIST(3)=LASTK
c      LAST=LAST+3
c      LAST=3
c      IF (LASTK .EQ. NF) THEN
c            LAST=LAST+1
c            MODHIST(LAST)=8HCOMPLETE
c            ENDIF
c      MODHIST(1)=LAST
      open (7,file='MODHIST',status='old')
      do
        read (7,'(A80)',end=11) card
      enddo
 11   continue
c      close (7)
C***  UPDATING THE MODEL HISTORY
      write(LCARD,88) JOBNUM,'. WRCONT   LASTK=',LASTK,TOC(TIMER),' sec'
   88 FORMAT (1H/,I3,A,I5,' COMPLETE  -  run time all K: ',i10,A)
      write (7,'(A100)') LCARD
      close (7)

c      IF (LASTK .EQ. NF) THEN
c            NEXTK=1
c         ELSE
c            NEXTK=LASTK+1
c         ENDIF
c      CALL WRITMS (3,NEXTK,1,5HNEXTK,-1,IDUMMY,IERR)

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
      ! CALL WRITMS (3,GTOT,ND,4HGTOT,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,ETOT,ND,4HETOT,-1,IDUMMY,IERR)
      ! CALL WRITMS (3,XTOT,ND,4HXTOT,-1,IDUMMY,IERR)
      !*** if POPNUM_CP is set, then check if writeout
      if(DECSTAR_OPT%POPNUM_CP>0) then
        if(mod(JOBNUM,DECSTAR_OPT%POPNUM_CP)==0) call cpPOPNUM(JOBNUM)
      endif
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
      ! CALL CLOSMS (3)
      ! CALL CLOSMS (7)
      ! CALL JSYMSET (2LG0,0)
!      close (6)
      ! STOP 'O.K.'
      RETURN
      END subroutine
      end module
