      MODULE MOD_ETL

      CONTAINS

      SUBROUTINE ETL (job)

      use UTILS, only: ASSERT
      use MOD_READMOD
      use MOD_WRITMS
      use MOD_WRITMSI
      use MOD_READPOP
      use MOD_COOP_M
      use MOD_READMS
      use MOD_READMSI
      use MOD_LIOP
      use MOD_DATOM_M
      use MOD_READRAD
      use MOD_ETLRAY
      use MOD_extUray
      use MOD_DIFFUS
      use MOD_ELIMIN
      use MOD_PRIINTL
      use MOD_STORXJL
      use MOD_FLGRID
      use MOD_DECETL
      use MOD_PRELINE
      use MOD_ETLHIST
      use MOD_WRITMOD
      use MOD_REBLANK
      use MOD_PRIETL
      use MOD_WRITRADL
      use MOD_FORMATS
      use MOD_ERROR
      use MOD_TICTOC
      use PARAMS_ARRAY
      use ABUNDANCES
      USE CONSTANTS

      USE COMMON_BLOCK
      USE FILE_OPERATIONS

      use VARDATOM

      IMPLICIT REAL*8(A - H, O - Z)

      INTEGER :: LASTUPD
      integer :: timer
      integer :: NDUMMY0
      
      real*8, allocatable :: DUMMY2(:,:)
      real*8 :: DUMMY0

!      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)

      COMMON /COMLBKG/ LBKG,XLBKG1,XLBKG2

      INTEGER XLBKG1,XLBKG2

      LOGICAL LBKG

!     RINAT TAGIROV:
!     See Fig. 7-29 in Mihalas, Stellar Atmospheres, 2nd edition, 1978
!     for the geometry of the ray-by-ray solution as well as
!     Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows.
!     I.Computational method for equivalent-two-level-atom source functions"
 
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
     $ ,EN(NDIMP2),V1(NDIMP2),V2(NDIMP2),RATCO(NDIMP2,NDIMP2)
     $ ,V4(NDIMP2),V5(NDIMP2),ENDELTA(NDIMP2),SCOLD(NFDIM,NDDIM)
     $ ,VOLD(NDIMP2),DB(NDIMP2,NDIMP2),WCHARM(NDDIM,NFDIM)
     $ ,DM(NDIMP2,NDIMP2),CRATE(NDIM,NDIM),RRATE(NDIM,NDIM)
     $ ,ETA(NDDIM),OPA(NDDIM),THOMSON(NDDIM),TAUTHOM(NDDIM)
     $ ,TNEW(NDDIM),XJCAPP(NFDIM),OPAC(NFDIM),SCNEW(NFDIM),DOPA(NFDIM)
     $ ,DETA(NFDIM),ETAC(NFDIM),EXPFAC(NFDIM),SIGMAKI(NFDIM,NDIM)
     $ ,PHI(NFLDIM),PWEIGHT(NFLDIM),DEPART(NDDIM,NDIM)
     $ ,OPAL(NDDIM),ETAL(NDDIM),ASF(NDDIM),BSF(NDDIM),PP(NDDIM)
     $ ,XJLMEAN(NDDIM),HBLUWI(NDDIM)
     $ ,TA(NDDIM),TB(NDDIM),TC(NDDIM),UB(NDDIM),GA(NDDIM),H(NDDIM)
     $ ,QQ(NDDIM),S(NDDIM),V(NDDIM),VA(NDDIM),VB(NDDIM)
     $ ,Uray(NDDIM)
     $ ,AF(NFLDIM,NFLDIM),BF(NFLDIM,NFLDIM)
     $ ,XJ(NFLDIM,NDDIM),XH(NFLDIM,NDDIM),XK(NFLDIM,NDDIM)
     $ ,XN(NFLDIM,NDDIM),W0(NDDIM),W1(NDDIM),W2(NDDIM),W3(NDDIM)
     $ ,BMHO(NFLDIM),BMNO(NFLDIM),BMHI(NFLDIM),BMNI(NFLDIM)
!     $ ,KODAT(MAXATOM),NFIRST(MAXATOM),NLAST(MAXATOM)
!     $ ,NCHARG(NDIM),MAINQN(NDIM),NOM(NDIM)
     $ ,ITNE(NDDIM),NFEDGE(NDIM)
     $ ,IWARN(NDDIM),LEVELPL(NDIM),MODHIST(MAXHIST)

      COMMON /COMELI/ XJCIND(NDDIM)

C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY !
     $ ,A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM),BX(NPDIM,NPDIM,NDDIM)
     $ ,WX(NPDIM,NDDIM)
!     $ ,INDNUP(MAXIND),INDLOW(MAXIND)

      LOGICAL :: PR_COND

      REAL*8, ALLOCATABLE :: U(:, :)

!     the Local approximate lambda-Operator for a given line
      REAL*8, ALLOCATABLE :: LO(:)

!      REAL*8, ALLOCATABLE :: AW(:, :)

!      REAL*8, ALLOCATABLE :: OpticalDepthLine(:)
!      REAL*8, ALLOCATABLE :: OpticalDepthLineCont(:)

      parameter (IPDIM = 25, NBDIM = 99)

      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)

      CHARACTER MODHEAD*104,MODHOLD*104, CARD*80, LCARD*420
!      CHARACTER LEVEL(NDIM)*10, CNAME*10, CREAD*10
      CHARACTER CNAME*10, CREAD*10
      CHARACTER*10 MAINPRO(NDDIM),MAINLEV(NDDIM)
!      CHARACTER*10 ELEMENT(MAXATOM), NAMEU, NAMEJ
      CHARACTER*10 NAMEU, NAMEJ
!      CHARACTER*4 KEYCOL(NDIM,NDIM)
!      CHARACTER*2 SYMBOL(MAXATOM)
      CHARACTER*7 JOB
      character*7 LINE(MAXIND)
      LOGICAL ETLKEY,LINEKEY(maxind)
      integer,external :: time
!      character*8 :: agaunt
      logical :: ierr9

      print *,'etl: Enter etl: '//writeTOC()
      call tic(timer)

C***  DECODING INPUT CARDS
      CALL DECETL(LSOPA,LSINT,VDOP,LINE,NLINE,LINEKEY,MAXIND,LBLANK)

      IF (NLINE .EQ. 0) THEN; WRITE(6, *) 'NO LINE OPTIONS DECODED'; GOTO 99; ENDIF
 
      CALL FLGRID(NFLDIM, NFL, PHI, PWEIGHT, DELTAX)

C***  changes by Margit Haberreiter
	CALL DATOM_M(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $               EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $               INDNUP,INDLOW,LASTIND,NATOM,
     $               ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $               NLAST,WAVARR,SIGARR,NFDIM)
	
C***  READING OF THE MODEL FILE
      IFL = 3

      open(IFL, file='MODFILE', STATUS='OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOPOLD,NF,XLAMBDA,FWEIGHT,AKEY,
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,
     $             NDDIM,NPDIM,NFDIM,LBLANK)

      close(IFL)

      print *,'ETL: ND = ',ND,' NP = ',NP
      call assert(ND*NP>1,'ND*NP <= 1')

!     LOO is declared in common_block.for
      IF (ALLOCATED(LOO)) DEALLOCATE(LOO)

      ALLOCATE(U(ND, NP))
      ALLOCATE(LOO(ND, LASTIND))

      IF (JOBNUM .GE. 1000) JOBNUM = JOBNUM - 100
      PRINT *,'ETL: VDOP=',VDOP
      If (vdop.le.0.) vdop=vdopold

      IFL=3
      open (IFL,file='POPNUM',STATUS='OLD')

      call readpop (ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum)

      close (ifl)

c***  read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             NCHARG,EDDARR,EDDI,NOM,WCHARM,N,lastind,
     $             EINST,MODHEAD,JOBNUM)

      JOBNUM = JOBNUM + 1

      if (vdop.ne.vdopold) then
         IFL=3 
         open (IFL,file='MODFILE',STATUS='OLD')
         print *,' Line Doppler width changed'
         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)
         close(IFL)
      endif
C***  NO LINE BLANKETING WITHIN A LINE 
c     switch off LB       
      LBLANK=0

      CALL REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

!     INTRODUCING DIMENSIONLESS VELOCITY UNITS:
!     ********************************************
!     RINAT TAGIROV:
!     The description of VDOP is given in LIOP.for
!     ********************************************
      DO L = 1, ND

         VELO(L) =  VELO(L)  / VDOP

         GRADI(L) = GRADI(L) / VDOP

      ENDDO
 
C***  CHECK WHETHER THESE DATA EXIST AND BELONG TO THE PRESENT MODEL
      IFL=9
	ierr9=.false.
      open (IFL,file='RADIOCL',STATUS='UNKNOWN')
      IERR=-10
      CNAME='MODHEAD'
      READ (ifl,'(A10)',err=1001,end=101) cread
      if (cread.eq.cname) then
         READ (ifl,'(A104)',err=1001,end=101) MODHOLD
         IF (MODHEAD .NE. MODHOLD) then
            goto 101
         else
            ierr=0
	   endif
	endif
101    continue

      IF (IERR.ne.0 .OR. JOB.EQ.'newline')THEN
         rewind (ifl)
         CNAME='MODHEAD'
         write (ifl,'(A10)')  cname
	   write (ifl,'(A104)') MODHEAD 
         CNAME='JOBNUM'
         CALL WRITMSI1(IFL,JOBNUM,CNAME,-1,IERR)
      ELSE
	   cname='JOBNUM'
         CALL READMSI1(ifl,LASTUPD,cname,ierr)
         IF (LASTUPD.GE. JOBNUM) PRINT *,
     $     ' ETL: WARNING - BACKGROUND DATA YOUNGER THAN PRESENT MODEL'
         IF (LASTUPD .LT. JOBNUM-20) PRINT *,
     $     ' ETL: WARNING - BACKGROUND DATA OLDER THAN 20 JOBS'
	   ierr9=.true.
      ENDIF
 
C***  NEWBGC = COUNTER OF LINES FOR WHICH NEW BACKGROUND CONT. RAD. FIELD IS
C***  CALCULATED
      NEWBGC = 0
 
C***  BEGINNING OF LOOP FOR EACH LINE ******************************************
C***  ONLY LINES WHICH APPEAR IN AN INPUT OPTION ARE TREATED !
C***  ( IN THE SEQUENCE OF THEIR OCCURENCE )

      ALLOCATE(LO(ND))

!      ALLOCATE(AW(NP, ND))

!      ALLOCATE(OpticalDepthLine(ND))
!      ALLOCATE(OpticalDepthLineCont(ND))


!      OPEN(UNIT = 14, FILE = 'line_feautrier_matrix.out',     FORM = 'FORMATTED')
!      OPEN(UNIT = 16, FILE = 'line_lambda_operator_diag.out', FORM = 'FORMATTED')

!      WRITE(14, 10231)
!      WRITE(16, 10233)

      SIGMAKI(1, 1) = -1.

      LineNumber = 1

      DO 7 NL = 1, NLINE

      ETLKEY = LINEKEY(NL)
 
C***  PREPARING SOME QUANTITIES FOR THE CONSIDERED LINE

      CALL PRELINE(NUP,LOW,IND,N,LRUD,XLAM,ND,NFL,LINE,BMHO,BMNO,
     $             BMHI,BMNI,XJLMEAN,HBLUWI,XJ,XH,XK,XN,ELEVEL,NL,
     $             NDIM,EINST,INDNUP,INDLOW,LASTIND)

      IF ( (NUP .EQ. 0) .OR. (LRUD .EQ. 0) ) GOTO 7 ! THIS IS HOW WE GET 75 LINES OUT OF NLINE = 150
 
C***  COMPUTATION OF THE ETLA SOURCE FUNCTION COEFFICIENTS, ASF AND BSF
C***  RATES AND POPNUMBERS ARE UPDATED !
C***  THIS ROUTINE IS SKIPPED FOR SUBSEQUENT "FORMAL" LINES
C***  EXCEPT OF THE FIRST FORMAL LINE WHICH FOLLOWS AN ETL LINE
      IF (ETLKEY) GOTO 2
      IF (NL .EQ. 1) GOTO 3
      IF (LINEKEY(NL-1)) GOTO 2
      GOTO 3
    2 CONTINUE
c***  ETLA option not active
         print *,' ETLA option is not active'
         write (6,*) ' ETLA option is not active'
	   stop 'ETLA'
    3 CONTINUE
 
C***  INNER BOUNDARY VALUES
      CALL DIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR)
 
C***  CONTINUUM AND LINE OPACITIES
      CALL COOP_M(XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $            OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $            N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $            ALPHA,SEXPO,AGAUNT,0,DUMMY2,WAVARR,SIGARR,
     $            LBKG,XLBKG1,XLBKG2,NF)

C***  BACKGROUND CONTINUUM RADIATION FIELD

      write (NAMEU,FMT_KEY) 'U   ',IND
      write (NAMEJ,FMT_KEY) 'XJC ',IND
      if (ierr9) then

        !*** READ U AND XJC FROM FILE 9
      IERR=1
      CALL READMS (ifl,U,ND*NP,NAMEU,IERR)
      call assert(ierr==0,'Error reading U from File')
      CALL READMS (ifl,XJCIND,ND,NAMEJ,ierr)
      call assert(ierr==0,'Error reading XJCIND from File')

      else

c234567890 234567890 234567890 234567890 234567890 234567890 234567890 2
C***     DATA FOR THAT LINE NOT PRESENT - CALCULATE IT NEW ...

      CALL ELIMIN(XLAM,DUMMY1,DUMMY0,U,Z,A,B,C,W,BX,WX,XJCIND,RADIUS,P,
     $            BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

c***     ... and write it for use in the next iteration
      CALL WRITMS (ifl,U,ND*NP,NAMEU,-1,IERR)
      CALL WRITMS (ifl,XJCIND,ND,NAMEJ,-1,IERR)

      NEWBGC = NEWBGC + 1

      ENDIF
 
      !***  ADD THE THOMSON EMISSIVITY, USING THE CONTINUUM RADIATION FIELD :
      !DO 9 L=1,ND
      ETA(:ND)=ETA(:ND)+OPA(:ND)*THOMSON(:ND)*XJCIND(:ND)
      !9 CONTINUE

      CALL LIOP_RTE(EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     $              ND,XLAM,ENTOT,POPNUM,RSTAR,OPAL,ETAL,VDOP, N)
 
C***  FORMAL SOLUTION OF RAD.TRANSFER IN THE COMOVING FRAME
C***  IN ORDER TO OBTAIN THE 0. TO 3. MOMENTS OF THE RADIATION FIELD
C***  AND HENCE THE EDDINGTON FACTORS

!      OpticalDepthLine(1) = 0.0D0
!      OpticalDepthLineCont(1) = 0.0D0

!      DO L = 2, ND

!         OpticalDepthLine(L) = OpticalDepthLine(L - 1) + 
!     $ ((OPAL(L) + OPAL(L - 1)) / 2.) * 
!     $ ((HEIGHT(L - 1) - HEIGHT(L)) / 
!     $ SolarRadiusKM)

!         OpticalDepthLineCont(L) = OpticalDepthLineCont(L - 1) +
!     $ ((OPA(L) + OPAL(L) + OPA(L - 1) + OPAL(L - 1)) / 2.0D0) *
!     $ ((HEIGHT(L - 1) - HEIGHT(L)) / 
!     $ SolarRadiusKM)

!      ENDDO

!      OpticalDepthLine(1) = OpticalDepthLine(2)
!      OpticalDepthLineCont(1) = OpticalDepthLineCont(2)

      LO(1 : ND) = 0.0D0

!      AW(1 : NP, 1 : ND) = 0.0D0

C***  LOOP FOR RAY-BY-RAY COMPUTATION

      DO 13 JP = 1, NP - 1

      LMAX = MIN0(NP + 1 - JP, ND)

      CALL EXTURAY(U, URAY, ND, NP, JP)

      CALL ETLRAY(URAY,Z,OPA,OPAL,ETA,ETAL,XJLMEAN,RADIUS,ND,NP,JP,P,
     $            UB,GA,H,QQ,S,V,VA,VB,VELO,GRADI,PP,
     $            BMHO,BMNO,BMHI,BMNI,BCORE,DBDR,
     $            XJCIND,ELEVEL,LOW,
     $            PHI,PWEIGHT,NFL,DELTAX,XJ,XH,XK,XN,W0,W1,W2,W3,
     $            NL, XLAM, LO)

!      AW(JP, 1 : LMAX) = W0(1 : LMAX)

   13 CONTINUE

!      do iii = 1, ND

!         print*, 'etl 0:', NL, iii, LO(iii)

!      enddo

!      stop

      LSOPAP = 1

      IF (LSOPA.GT.0 .and. nl.eq.LSOPA)
     $      CALL PRIETL (IND,XLAM,ND,OPA,OPAL,ETA,ETAL,RADIUS,
     $      JOBNUM,LSOPAP,XJLMEAN,MODHEAD,ASF,BSF,T,ETLKEY)
      IF (.NOT.ETLKEY) GOTO 20
         print *,' ETLA option is not active'
         write (6,*) ' ETLA option is not active'
	   stop 'ETLA'
 
C***  UPDATING THE LINE RADIATION FIELD IN THE MODEL FILE AND IN THE CURRENT XJL
   20 continue

!     RINAT TAGIROV:
!     We found an abnormal behavior of the local operator element
!     at the boundary points: LO(ND - 1) > LO(ND) and LO(1) > LO(2).
!     It might have something to do with the boundary conditions in the
!     Feautrier system for the solution of the radiative transfer equation (see CMFSET.FOR).
!     See Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows.
!     I.Computational method for equivalent-two-level-atom source functions" for boundary conditions formulation.
!     It follows from the paper that the boundary conditions are of the first order only.
!     It is consistent with the comments in CMFSET.FOR (see the core rays for the boundary
!     condition at the innermost point of the atmosphere model)
!     This might be the cause of such behavior.
!     In the paper it says that refining the atmosphere model
!     near the boundary points can help. I did that but no improvement was attained.
!     Maybe I did the refinement of the atmosphere model wrongly.
!     Another idea that comes to mind is that this behavior after all may not have
!     anything to do with the boundary conditions, but with the angular and
!     frequency weights W0 and PWEIGHT (see ETLRAY.FOR).
!     If there is something wrong with them it might affect both intensity XJLMEAN and LO element,
!     but I did not check this possibility.
!     Be that as it may, in what follows, to obviate this irregular behavior I linearly
!     extrapolated LO(ND) and LO(1) using the two preceding LO elements.

      LO(1) =  EXTRAP_LO_B_VAL(LO, ND, 3, 2, 1)
      LO(ND) = EXTRAP_LO_B_VAL(LO, ND, ND - 2, ND - 1, ND)

!     XJLMEAN is stored in XJL, LO is stored in LOO (which means LO OVERALL, i.e. for all lines at once)
      CALL STORXJL(XJL, XJLMEAN, ND, LASTIND, IND, LOO, LO)

      LineNumber = LineNumber + 1

    7 CONTINUE
!     END OF LOOP FOR EACH LINE

!      do iii = 1, LASTIND

!         print*, 'etl:', iii, LOO(L, iii)

!      enddo

!      stop

      DEALLOCATE(LO)

!      DEALLOCATE(AW)

!      DEALLOCATE(OpticalDepthLine)
!      DEALLOCATE(OpticalDepthLineCont)

      CLOSE(ifl)

!      CLOSE(14)

c***  store the line radiation field in file RADIOL
      call writradl(XJL,XJLMEAN,EINST,NCHARG,NOM,ND,N,LASTIND,MODHEAD,JOBNUM)

      IF (LSINT.GT.0) CALL PRIINTL(NDIM,N,LEVEL,WEIGHT,EINST,LASTIND,LINE,NLINE,
     $                             INDLOW,INDNUP,ELEVEL,ND,XJL,LSINT,JOBNUM,MODHEAD)

C***  UPDATING THE MODEL HISTORY
      CALL ETLHIST(NLINE,LINE,LINEKEY,MODHIST,MAXHIST,JOBNUM,VDOP,VDOPOLD,LAST,LCARD)

      open (7,file='MODHIST',status='old')
  8   read (7,'(A80)',end=11) card
	goto 8
 11   continue
	if (last.gt.400) last=400
      write (LCARD(LAST+1:LAST+20),'(" Rtime:",i9," sec")') toc()
      write (7,'(A120)') LCARD
      close (7)
 
C***  LOGFILE REMARK ON CONTINUUM BACKGROUND CALCULATIONS
      IF (NEWBGC .GT. 0) THEN
	   print '("NEW BACKGROUND CONTINUUM CALCULATED FOR"'//
     $          ',I4," LINES")',newBGC
      ENDIF
 
   99 CONTINUE

      close (6)
      RETURN

1001  continue
      write (6,*) ' ETL: error exit 1001'
	call error('etl: error exit 1001')

10231 FORMAT(2x,'NL',7x,'JP',9x,'K',8x,'ID',14x,
     $ 'TAL',20x,'TBL',20x,'TCL',16x,'LambdaOperDiag')

!10233 FORMAT('LN',14x,'IND',13x,'XLAM',13x,'ID',8x,'HEIGHT',8x,
!     $ 'OpacLine',8x,'EmisLine',8x,'OptDepLine',8x,'LamOpDiagMean')

!10234 FORMAT(I2,5x,I3,5x,E15.7,5x,I2,5x,E15.7,5x,E15.7,5x,E15.7,5x,E15.7,5x,E23.15)

      END SUBROUTINE


      FUNCTION EXTRAP_LO_B_VAL(LO, LO_SIZE, PD, UD, BD) RESULT(LO_B_VAL)

!     LINEARLY EXTRAPOLATES THE LOCAL OPERATOR ELEMENT BOUNDARY VALUE FROM THE TWO PRECEDING VALUES

      USE COMMON_BLOCK

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::                    LO_SIZE

      REAL*8, DIMENSION(LO_SIZE), INTENT(IN) :: LO

      INTEGER, INTENT(IN) ::                    BD, UD, PD

      REAL*8 ::                                 C1, C2

      REAL*8 ::                                 LO_B_VAL

      C1 = (LO(PD) - LO(UD)) / (HEIGHT(PD) - HEIGHT(UD))

      C2 = LO(UD) - C1 * HEIGHT(UD)

      LO_B_VAL = C1 * HEIGHT(BD) + C2

      RETURN

      END FUNCTION EXTRAP_LO_B_VAL

      END MODULE
