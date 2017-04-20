      MODULE MOD_ETL

      CONTAINS

      SUBROUTINE ETL(job)

      use UTILS, only: ASSERT
      use MOD_READMOD
      use MOD_WRITMS
      use MOD_WRITMSI
      use MOD_READPOP
      use MOD_COOP
      use MOD_READMS
      use MOD_READMSI
      use MOD_LIOP
      use MOD_DATOM
      use MOD_READRAD
      use MOD_ETLRAY
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
      USE CONSTANTS

      use vardatom_lte
      use vardatom_nlte
      use varhminus
      use varsteal
      use common_block
      use file_operations
      use math
      use phys
      use local_operator

      implicit real*8(a - h, o - z)

      INTEGER :: LASTUPD
      integer :: timer
      integer :: NDUMMY0
      
      real*8, allocatable :: DUMMY2(:, :)
      real*8 :: DUMMY0

      COMMON /COMLBKG/ LBKG, XLBKG1, XLBKG2

      INTEGER XLBKG1, XLBKG2

      LOGICAL LBKG

!     RINAT TAGIROV:
!     See Fig. 7-29 in Mihalas, Stellar Atmospheres, 2nd edition, 1978
!     for the geometry of the ray-by-ray solution as well as
!     Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows.
!     I.Computational method for equivalent-two-level-atom source functions"

!     the Local approximate lambda-Operator for a given line
      REAL*8, ALLOCATABLE :: LO(:)

      CHARACTER MODHEAD*104,MODHOLD*104, CARD*80, LCARD*420
      CHARACTER CNAME*10, CREAD*10
      CHARACTER*10 NAMEU, NAMEJ
      CHARACTER*7 JOB
      character*7 LINE(lastind_nlte)
      LOGICAL ETLKEY, LINEKEY(lastind_nlte)
      logical :: ierr9

      real*8, allocatable, dimension(:) ::    opal, etal

      real*8, allocatable, dimension(:) ::    xjcind, xjlmean

      print*, 'entered etl: ' // writeTOC()
      call tic(timer)

!     DECODING INPUT CARDS
      CALL DECETL(LSOPA,LSINT,VDOP,LINE,NLINE,LINEKEY,lastind_nlte,LBLANK)

      IF (NLINE .EQ. 0) THEN; WRITE(6, *) 'NO LINE OPTIONS DECODED'; GOTO 99; ENDIF
 
      CALL FLGRID(NFL, PHI, PWEIGHT, DELTAX)

      CALL DATOM(datom_lte,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $           EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $           INDNUP,INDLOW,LASTIND,NATOM,
     $           ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $           NLAST,WAVARR,SIGARR,NFDIM)
	
!***  READING OF THE MODEL FILE
      IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'OLD')

      CALL READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $             GRADI,RSTAR,VDOPOLD,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      call assert(ND > 1, 'ND <= 1'); call assert(NP > 1, 'NP <= 1')

      if (allocated(llo))       deallocate(llo);       allocate(llo(ND, lastind_nlte))
      if (allocated(tau_line))  deallocate(tau_line);  allocate(tau_line(ND, lastind_nlte))
      if (allocated(damp_line)) deallocate(damp_line); allocate(damp_line(ND, lastind_nlte))

      if (allocated(opal))      deallocate(opal);      allocate(opal(ND))
      if (allocated(etal))      deallocate(etal);      allocate(etal(ND))

      if (allocated(xjcind))    deallocate(xjcind);    allocate(xjcind(ND))
      if (allocated(xjlmean))   deallocate(xjlmean);   allocate(xjlmean(ND))

      damp_line(1 : ND, 1 : lastind_nlte) = .false.

      llo(1 : ND, 1 : lastind_nlte) = 0.0d0

      If (vdop .le. 0.) vdop = vdopold

      IFL = 3; open(IFL, file = 'POPNUM', STATUS = 'OLD')

      call readpop(ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum)

      close (ifl)

!***  read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $             ncharg_nlte,EDDARR,EDDI,nom_nlte,WCHARM,N_nlte,lastind_nlte,
     $             einst_nlte,MODHEAD,JOBNUM)

      JOBNUM = JOBNUM + 1

      if (vdop.ne.vdopold) then
         IFL=3 
         open (IFL,file='MODFILE',STATUS='OLD')
         print *,' Line Doppler width changed'
         CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                GRADI,RSTAR,VDOP,NF,
     $                XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $                ABXYZ,NATOM,MODHEAD,JOBNUM)
         close(IFL)
      endif
!***  NO LINE BLANKETING WITHIN A LINE 
!     switch off LB       
      LBLANK=0

      CALL REBLANK(LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

!     INTRODUCING DIMENSIONLESS VELOCITY UNITS:
!     ********************************************
!     RINAT TAGIROV:
!     The description of VDOP is given in LIOP.for
!     ********************************************
      DO L = 1, ND

         VELO(L) =  VELO(L)  / VDOP

         GRADI(L) = GRADI(L) / VDOP

      ENDDO
 
!***  CHECK WHETHER THESE DATA EXIST AND BELONG TO THE PRESENT MODEL
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
 
!***  NEWBGC = COUNTER OF LINES FOR WHICH NEW BACKGROUND CONT. RAD. FIELD IS
!***  CALCULATED
      NEWBGC = 0
 
!***  BEGINNING OF LOOP FOR EACH LINE ******************************************
!***  ONLY LINES WHICH APPEAR IN AN INPUT OPTION ARE TREATED !
!***  ( IN THE SEQUENCE OF THEIR OCCURENCE )

      ALLOCATE(LO(ND))

      SIGMAKI(1, 1) = -1.

      LineNumber = 1

      DO 7 NL = 1, NLINE

!      lte_line = .false.

      ETLKEY = LINEKEY(NL)
 
!***  PREPARING SOME QUANTITIES FOR THE CONSIDERED LINE

      CALL PRELINE(NUP,LOW,IND,N_nlte,LRUD,XLAM,ND,NFL,LINE,
     $             XJLMEAN,elevel_nlte,NL,einst_nlte,
     $             indnup_nlte,indlow_nlte,lastind_nlte)

!      IF ((NUP .EQ. 0) .OR. (LRUD .EQ. 0) .or. lte_line) GOTO 7
      IF ((NUP .EQ. 0) .OR. (LRUD .EQ. 0)) GOTO 7

!      print*, 'etl: NL = ', NL, ' out of ', NLINE

      print*, 'etl: LN = ', LineNumber
 
!***  RATES AND POPNUMBERS ARE UPDATED !
!***  THIS ROUTINE IS SKIPPED FOR SUBSEQUENT "FORMAL" LINES
!***  EXCEPT OF THE FIRST FORMAL LINE WHICH FOLLOWS AN ETL LINE
      IF (ETLKEY) GOTO 2
      IF (NL .EQ. 1) GOTO 3
      IF (LINEKEY(NL-1)) GOTO 2
      GOTO 3
    2 CONTINUE
!***  ETLA option not active
         print *,' ETLA option is not active'
         write (6,*) ' ETLA option is not active'
	   stop 'ETLA'
    3 CONTINUE
 
!***  INNER BOUNDARY VALUES
      CALL DIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR)
 
!***  CONTINUUM AND LINE OPACITIES
      CALL COOP(XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $          OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $          N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $          ALPHA,SEXPO,AGAUNT,0,DUMMY2,
     $          WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF),
     $          LBKG,XLBKG1,XLBKG2,NF)

!***  BACKGROUND CONTINUUM RADIATION FIELD

      write (NAMEU,FMT_KEY) 'U   ',IND
      write (NAMEJ,FMT_KEY) 'XJC ',IND
      if (ierr9) then

!     READ U AND XJC FROM FILE 9
      IERR=1
      CALL READMS (ifl,U,ND*NP,NAMEU,IERR)
      call assert(ierr==0,'Error reading U from File')
      CALL READMS (ifl,XJCIND,ND,NAMEJ,ierr)
      call assert(ierr==0,'Error reading XJCIND from File')

      else

!***     DATA FOR THAT LINE NOT PRESENT - CALCULATE IT NEW ...

      CALL ELIMIN(XLAM,DUMMY1,DUMMY0,U,Z,XJCIND,RADIUS,P,
     $            BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

!***     ... and write it for use in the next iteration
      CALL WRITMS(ifl,U,ND*NP,NAMEU,-1,IERR)
      CALL WRITMS(ifl,XJCIND,ND,NAMEJ,-1,IERR)

      NEWBGC = NEWBGC + 1

      ENDIF
 
      !ADD THE THOMSON EMISSIVITY, USING THE CONTINUUM RADIATION FIELD :
      ETA(:ND)=ETA(:ND)+OPA(:ND)*THOMSON(:ND)*XJCIND(:ND)

      CALL LIOP_RTE(EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     $              ND,XLAM,ENTOT,POPNUM,RSTAR,OPAL,ETAL,VDOP, N)

      tau_line(1 : ND, ind) = opt_dep(opal(1 : ND) / RSTAR, 1.0D+5 * height(1 : ND), ND)
 
!***  FORMAL SOLUTION OF RAD.TRANSFER IN THE COMOVING FRAME
!***  IN ORDER TO OBTAIN THE 0. TO 3. MOMENTS OF THE RADIATION FIELD
!***  AND HENCE THE EDDINGTON FACTORS

      LO(1 : ND) = 0.0D0

!***  LOOP FOR RAY-BY-RAY COMPUTATION

      DO 13 JP = 1, NP - 1

      LMAX = MIN0(NP + 1 - JP, ND)

      CALL ETLRAY(U(1 : ND, JP),Z,OPA,OPAL,ETA,ETAL,XJLMEAN,RADIUS,ND,NP,JP,P,
     $            VELO,GRADI,BCORE,DBDR,PHI,PWEIGHT,NFL,DELTAX,LO)

   13 CONTINUE

      LSOPAP = 1

      IF (LSOPA.GT.0 .and. nl.eq.LSOPA)
     $      CALL PRIETL(IND,XLAM,ND,OPA,OPAL,ETA,ETAL,RADIUS,
     $      JOBNUM,LSOPAP,XJLMEAN,MODHEAD,T,ETLKEY)
      IF (.NOT.ETLKEY) GOTO 20
         print *,' ETLA option is not active'
         write (6,*) ' ETLA option is not active'
	   stop 'ETLA'
 
!***  UPDATING THE LINE RADIATION FIELD IN THE MODEL FILE AND IN THE CURRENT XJL
   20 continue

!     RINAT TAGIROV:
!     We found an abnormal behavior of the local operator element
!     at the boundary points: LO(ND - 1) > LO(ND) and LO(1) > LO(2).
!     It might have something to do with the boundary conditions in the
!     Feautrier system for the solution of the radiative transfer equation (see CMFSET).
!     See Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows.
!     I.Computational method for equivalent-two-level-atom source functions" for boundary conditions formulation.
!     It follows from the paper that the boundary conditions are of the first order only.
!     It is consistent with the comments in CMFSET (see the core rays for the boundary
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

      lo(1) =  extrap_to_boundary(ND, height, lo, 3,      2,      1)
      lo(ND) = extrap_to_boundary(ND, height, lo, ND - 2, ND - 1, ND)

!     XJLMEAN is stored in XJL, LO is stored for each line in LLO (which means Line Local Operator)
      CALL STORXJL(XJL, XJLMEAN, ND, lastind_nlte, IND, llo, lo)

      LineNumber = LineNumber + 1

    7 CONTINUE
!     END OF LOOP FOR EACH LINE

!      stop

      deallocate(lo)
      deallocate(xjcind)

      close(ifl)

!     perform the acceleration damping in case the density is too high at the outer edge of the atmosphere
      if (damp_acc) call acc_damp(ND, lastind_nlte, tau_line, LLO, damp_line)

!***  store the line radiation field in file RADIOL
      call writradl(XJL,XJLMEAN,einst_nlte,ncharg_nlte,nom_nlte,ND,N_nlte,lastind_nlte,MODHEAD,JOBNUM)

      deallocate(xjlmean)

      IF (LSINT.GT.0) CALL PRIINTL(N_nlte,level_nlte,weight_nlte,einst_nlte,lastind_nlte,LINE,NLINE,
     $                             indlow_nlte,indnup_nlte,elevel_nlte,ND,XJL,LSINT,JOBNUM,MODHEAD)

!***  UPDATING THE MODEL HISTORY
      CALL ETLHIST(NLINE,LINE,LINEKEY,JOBNUM,VDOP,VDOPOLD,LAST,LCARD)

      open (7,file='MODHIST',status='old')
  8   read (7,'(A80)',end=11) card
	goto 8
 11   continue
	if (last.gt.400) last=400
      write (LCARD(LAST+1:LAST+20),'(" Rtime:",i9," sec")') toc()
      write (7,'(A120)') LCARD
      close (7)
 
!***  LOGFILE REMARK ON CONTINUUM BACKGROUND CALCULATIONS
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

10231 FORMAT(2x,'NL',7x,'JP',9x,'K',8x,'ID',14x, 'TAL',20x,'TBL',20x,'TCL',16x,'LambdaOperDiag')

      END SUBROUTINE

      END MODULE
