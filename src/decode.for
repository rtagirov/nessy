      module mod_decode

      type :: DECSTAR_PARAMS 
      integer :: POPNUM_CP
      end type
      type(DECSTAR_PARAMS) :: DECSTAR_OPT = DECSTAR_PARAMS(0)

      contains

      subroutine decstar(MODHEAD,FM,RSTAR,t_eff,glog,xmass,VDOP,TTABLE,TPLOT,NATOM,KODAT,IDAT,LBLANK,ATMEAN)

!     DECODES INPUT CARDS, CALLED FROM WRSTART AND WRCONT

      use common_block
      use file_operations
 
      IMPLICIT REAL*8(A - H, O - Z)

      PARAMETER (ONE = 1.D+0)
      
C***  COMMON/VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS TO FUNCTION WRVEL
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE
C***  COMMON /COMRPAR/ TRANSFERS RADIUS-GRID OPTION CARD TO SUBR. RGRID
      COMMON /COMRPAR/ RPAR
C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE TO SUBR. GREY AND PRIMOD
C***  ATTENTION: TRANSFER OF TEFF ALSO TO SUBR. DECFREQ, DTDR, JSTART !!!!!!!!!
C***             AND TO MAIN PROGRAM WRSTART !!!!!!!!!!!!!!!
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  COMMON / COMPLOT / TRANSFERS THE MAXIMUM OF THE Y-AXIS TO SUBR. PLOTT
      COMMON / COMPLOT / YMAX
      DIMENSION KODAT(NATOM)
!      dimension ABXYZ(NATOM)
      LOGICAL TTABLE, SPHERIC, TPLOT, FAL

      CHARACTER MODHEAD*104, KARTE*80,DAT*8, TIM*8
      CHARACTER RPAR*80
      REAL*8 ATMEAN
      
      TYPE ABUNDANCE_
                CHARACTER(4) name
                 INTEGER index
      END TYPE
       TYPE (ABUNDANCE_) :: ABUNDANCE

       DIMENSION ABUNDANCE(29)
! the abundances must be set here
       DATA ABUNDANCE / 
     &  Abundance_('ABH=',     2    ),
     &  Abundance_('ALI=',     3    ),
     &  Abundance_('ABE=',     4    ),
     &  Abundance_('ABB=',     5    ),
     &  Abundance_('ABC=',     6    ),
     &  Abundance_('ABN=',     7    ),
     &  Abundance_('ABO=',     8    ),
     &  Abundance_('ABF=',     9    ), 
     &  Abundance_('ANE=',    10    ), 
     &  Abundance_('ANA=',    11    ), 
     &  Abundance_('AMG=',    12    ),
     &  Abundance_('AAL=',    13    ),
     &  Abundance_('ASI=',    14    ),
     &  Abundance_('ABP=',    15    ),
     &  Abundance_('ABS=',    16    ),
     &  Abundance_('ACL=',    17    ),
     &  Abundance_('AAR=',    18    ),
     &  Abundance_('ABK=',    19    ),
     &  Abundance_('ACA=',    20    ),
     &  Abundance_('ASC=',    21    ),
     &  Abundance_('ATI=',    22    ),
     &  Abundance_('ABV=',    23    ),
     &  Abundance_('ACR=',    24    ),
     &  Abundance_('AMN=',    25    ),
     &  Abundance_('AFE=',    26    ), 
     &  Abundance_('ACO=',    27    ),
     &  Abundance_('ANI=',    28    ),
     &  Abundance_('ACU=',    29    ),
     &  Abundance_('AZN=',    30    ) /
     
C***  SOLAR RADIUS ( CM )
      DATA RSUN /6.96d10/

      XLBKG1 = 0
      XLBKG2 = 0

      RPAR=' '
      RSTAR=-ONE
      RMAX=-ONE
      VDOP=-ONE
      IDAT=-1
      LBLANK=0
      FM=.0
      XMDOT=.0
      ABXYZ(:)=0.
      FAL=.FALSE.
      TTABLE=.FALSE.
      SPHERIC=.FALSE.
      LBKG = .FALSE.
      TPLOT=.FALSE.
      YMAX=0.0
      TMODIFY=0.0
      TMIN=.0

C***  MODEL HEADER
      CALL DATE_AND_TIME (DAT)
      CALL CLOCK(TIM)
      MODHEAD='MODEL START'
      MODHEAD(15:)=DAT
      MODHEAD(25:)=TIM

      WRITE(*, '(/,A,/)') '------------------------------ CARDS FILE START ------------------------------'

!      OPEN (1, FILE = 'CARDS', STATUS = 'OLD'); REWIND 1
      OPEN (1, FILE = 'cards.inp', STATUS = 'OLD'); REWIND 1

   10 READ (1,11,END=66) KARTE
   11 FORMAT (A)
!      PRINT 2,KARTE
!    2 FORMAT (1X,A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
            ELSE
               PRINT *,' NO VALID OPTION'
               LBLANK=0
            ENDIF

            GOTO 10
            ENDIF

      IF (KARTE(:11) .EQ. 'PRINT DATOM' ) THEN

        IDAT=1
        GOTO 10
      ENDIF

!     FALP IS THE LOGICAL KEY FOR THE RISING TEMPERATURE PROFILE
      IF (KARTE(:3).EQ.'FAL') THEN

        FAL=.TRUE.
        GOTO 10
      ENDIF
      IF (KARTE(:4).EQ.'LBKG') THEN

! see comblock.for for the description of LBKG, XLBKG1, XLBKG2

        LBKG=.TRUE.

        DECODE (80,100,KARTE) XLBKG1,XLBKG2
  100   FORMAT (5X,I4,1X,I5)
        PRINT *,'LBKG= ',LBKG,XLBKG1,XLBKG2
        IF ((XLBKG1 .eq. 0) .or. (XLBKG1 .eq. 0) .or.
     $              (xlbkg2 .lt. xlbkg1)) THEN
          write (6,*) 'something wrong with lbkg-wavelength range'
          write (6,*) 'decstar: LBKG= ',LBKG,XLBKG1,XLBKG2
          print *, 'possibly wrong format lbkg-wavelength range'
          print *, 'decstar: LBKG= ',LBKG,XLBKG1,XLBKG2
          pause
        ENDIF
        GOTO 10
      ENDIF
      IF (KARTE(:5).EQ.'TABLE') THEN

            TTABLE=.TRUE.
            GOTO 10
            ENDIF
      IF (KARTE(:5).EQ.'SPHER') THEN

            SPHERIC=.TRUE.
            GOTO 10
            ENDIF
      IF (KARTE(:6).EQ.'PLOT T') THEN
C                       ======
            DECODE (80,6,KARTE) YMAX
    6       FORMAT (22X,F6.0)
            TPLOT=.TRUE.
            GOTO 10
            ENDIF
      IF (KARTE(:7).EQ.'TMODIFY') THEN
C                       =======
            DECODE (80,9,KARTE) TMODIFY
    9       FORMAT (8X,F10.0)
            GOTO 10
            ENDIF
      IF (KARTE(:4) .EQ. 'TMIN') THEN
C                         ====
            DECODE (80,8,KARTE) TMIN
    8       FORMAT (5X,F10.0)
            GOTO 10
            ENDIF
      IF (KARTE(:5).EQ.'RGRID') THEN
C                       =====
            RPAR=KARTE
            GOTO 10
            ENDIF
      IF (KARTE(:9).NE.'MASS FLUX') GOTO 14
C                       =========
      DECODE (80,13,KARTE) FM
   13 FORMAT (29X,F10.0)
      GOTO 10
     
   14 IF (KARTE(:5).EQ.'RSTAR') THEN
C                       =====
      DECODE (80,16,KARTE) RSTAR
   16 FORMAT (20X,F10.0)
      GOTO 10
      ENDIF


!================================= Abundancies ================================================

       DO I=1,SIZE(ABUNDANCE)

         IF (Karte(:4) .eq. (Abundance(I)%NAME)) THEN

           READ (KARTE,'(4X,G10.0)') ABNDNC
               

           IF ((KODAT(Abundance(I)%Index).LE.0.).AND.(ABNDNC.NE.0.))THEN
             PRINT *, 'INCONSISTENCY BETWEEN ATOMIC DATA AND ABUNDANCE'
             PRINT *,KARTE,ABNDNC
             PAUSE
             STOP 'ERROR'
           ELSE



             ! Set the given abundance at the given index
             ABXYZ(KODAT(Abundance(I)%Index))=ABNDNC

           ENDIF
           GOTO 10 ! go back to the beginning of the read file
         ENDIF
       END DO
      
! ===================== END ABUNDANCIES ================================

      IF (KARTE(:7).EQ.'ATMEAN=') THEN
C                       ====
            DECODE (80,'(8X,F10.0)',KARTE) ATMEAN
            
            ENDIF     


      IF (KARTE(:6).NE.'VELPAR') GOTO 20
C                       ======
      DECODE (80,19,KARTE) VFINAL,VMIN,BETA,RMAX
   19 FORMAT(22X,F7.1,5X,F6.0,5X,F4.0,5X,F6.0)
      GOTO 10
     
   20 IF (KARTE(:5).NE.'VDOP=') GOTO 23
C                       =====
      DECODE (80,22,KARTE) VDOP
   22 FORMAT (5X,F10.0)
      GOTO 10
     
   23 IF (KARTE(:8).NE.'HEADLINE') GOTO 27
C                        ==========
      MODHEAD(35:)=KARTE(10:)
      GOTO 10
     
   27 IF (KARTE(:5).NE.'MDOT=') GOTO 30
C                       =====
      DECODE (80,29,KARTE) XMDOT
   29 FORMAT (5X,F10.0)
      GOTO 10
     
   30 IF (KARTE(:5).NE.'TEFF=') GOTO 33
C                       =====
      DECODE (80,32,KARTE) TEFF
   32 FORMAT (5X,F10.0)
      GOTO 10
     
   33 IF (KARTE(:6).NE.'LOG G=') GOTO 36
C                       ======
      DECODE (80,35,KARTE) GLOG,XMASS
   35 FORMAT (6X,F10.0,6X,F10.0)
      GOTO 10
   36 IF(KARTE(:9).EQ.'POPNUM CP') then
        READ(KARTE,'(10XI8)') DECSTAR_OPT%POPNUM_CP
        PRINT *,'POPNUM CP = ',DECSTAR_OPT%POPNUM_CP
        GOTO 10
      ENDIF
      GOTO 10
     
   66 CONTINUE
      CLOSE (1)

      WRITE(*, '(/,A,/)') '------------------------------ CARDS FILE END ------------------------------'

C***  CHECK OF MISSING SPECIFICATIONS
      IF (FM.EQ..0 .AND. XMDOT .EQ. .0) THEN
            print*, 'MASS LOSS NOT SPECIFIED'
            PAUSE
            STOP "ERROR"
            ENDIF
      IF (RPAR .EQ. ' ') THEN
            print*, 'RADIUS GRID NOT SPECIFIED'
            PAUSE
            STOP "ERROR"
            ENDIF
      IF (RSTAR .LT. .0) THEN
            print*, 'RSTAR NOT SPECIFIED'
            PAUSE
            STOP "ERROR"
            ENDIF
      IF (RMAX .LT. .0) THEN
            print*, 'RMAX NOT SPECIFIED'
            PAUSE
            STOP "ERROR"
            ENDIF
      IF (VDOP .LT. .0) THEN
            print*,'VDOP IN CARDS NOT SPECIFIED'
            PAUSE "ERROR"
            STOP "ERROR"
            ENDIF
CMH    ATMEAN	 HAS TO BE AROUND 1
      IF (ATMEAN.LT. 1.e-1) THEN
            print*,'ATMEAN IN CARDS NOT DEFINED!'
            PAUSE "ERROR"
            STOP "ERROR"
            ENDIF
     
C***  COMPUTATION OF THE RELATIVE ABUNDANCE (BY NUMBER) OF HELIUM

      ABREST=0.
      DO 93 NA=2,NATOM
           IF (KODAT(NA) .GT. 0) then

           ABREST=ABREST+ABXYZ(KODAT(NA)) 

           endif
   93 CONTINUE

      IF (ABREST .GT. ONE) THEN
         WRITE (6,*) 'WRONG REL. ABUNDANCES', abrest
         PAUSE
         STOP 'ERROR'
         ENDIF
      IF (KODAT(1) .GT. 0) ABXYZ(KODAT(1))=ONE-ABREST
     
C***  CHECK OF UNNECESSARY ATOMIC DATA (WARNING!)
      DO 95 NA=1,NATOM
      IF (KODAT(NA) .GT. 0) THEN
         IF (ABXYZ(KODAT(NA)) .EQ. 0.) THEN
         PRINT *,' NA ',NA,KODAT(NA),ABXYZ(KODAT(NA))
         PRINT 96
   96    FORMAT (/////,10X,10('='),' WARNING: UNNECESSARY ATOMIC DATA ',
     $          'IN FILE DATOM DECODED ',10('='))
         WRITE (6,*) 'UNNECESSARY ATOMIC DATA DECODED'
         PAUSE
         STOP 'WARNING'
         ENDIF
      ENDIF
   95 CONTINUE
     
C***  COMPUTATION OF THE MASS FLUX
      IF (FM .NE. .0) FM= 10.**FM
c     3.01516 = log( Mo/(yr/sec)/Ro/Ro/4/Pi )
c     (this was just 3.02 !)
      IF (XMDOT .NE. .0) FM= 10.**(XMDOT+3.01516) /RSTAR/RSTAR

C***  STELLAR RADIUS IN CM
      RSTAR = RSTAR * RSUN

      t_eff = TEFF

      return

      end subroutine

      subroutine decste(LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,IFRRA,ITORA,LSEXPO,
     $                  IFLUX,IDAT,LEVELPL,N,IPLOTF,NEWWRC,
     $                  NGAMR,NGAML,AGAMR,AGAML,LINE,LASTIND,
     $                  TPLOT,Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB,LBLANK)

      use vardatom_nlte
      use common_block
      use file_operations

!     DECODING INPUT OPTIONS FOR STEAL.FOR

      IMPLICIT REAL*8(A - H, O - Z)

      CHARACTER KARTE*80
      DIMENSION LEVELPL(N)
      INTEGER   NGAMR(10),NGAML(10)
      DIMENSION AGAMR(10),AGAML(10)
      LOGICAL   LINE(LASTIND),PROLIB,NOTEMP,TPLOT
C***  COMMON / COMPLOT / TRANSFERS THE MAXIMUM OF THE Y-AXIS TO SUBR. PLOTT
      COMMON / COMPLOT / YMAX
     
C***  NO PRINT OPTIONS SET (DEFAULT)
      TEFFE=-1.
      LSRAT=-1
      LSPOP=-1
      LSEXPO=-1
      IHIST=-1
      IFLUX=-1
      IDAT=-1
      LBLANK=0
C***  OLDSTART=.F.:  POPNUMBERS ARE CALCULATED IN SUBR. POPZERO (DEFAULT)
C***  OLDSTART=.T.:  POPNUMBERS FROM OLD MODEL FILE
      OLDSTART=.FALSE.
C***  JOBMAX EXCEEDED (DEFAULT)
      JOBMAX=-1
C***  DEMANDED CONVERGENCE ACCURACY (DEFAULT)
      EPSILON=.005
C*** NO REDUCTION OF CORRECTIONS (DEFAULT)
      REDUCE = 1.
C***  NO POPNUMBER PLOT (DEFAULT)
      LEVELPL(1)=0
      NPLOT=0
C***  NO FLUX PLOT (DEFAULT)
      IPLOTF=-1
      PROLIB=.FALSE.
C***  NUMBER OF REPEAT ITERATIONS BETWEEN TWO WRCONT-JOBS (DEFAULT)
      NEWWRC=11
C***  NO SCHARMER AMPLIFICATION OF CORRECTIONS
      DO 21 I=1,10
      NGAMR(I)=999999999
      NGAML(I)=999999999
      AGAMR(I)=0.
      AGAML(I)=0.
   21 CONTINUE
      NGAMR(1)=0
      NGAML(1)=0
      MGC=1
      MGR=1
      MGL=1
      IPRICC=0
      IPRILC=0
C***  LINES THAT HAVE NOT BEEN TREATED IN THE RADIATION TRANSFER CALCULATIONS,
C***  I.E. WHICH ARE NOT QUOTED IN A 'LINE' OPTION CARD, ARE NOT SCHARMER-
C*** AMPLIFIED.
      DO 32 IN=1, LASTIND
      LINE(IN)=.FALSE.
   32 CONTINUE
C***  DISTANCE MODULUS (DEFAULT VALUE)
      Y0=.0
C***  NO SUPPRESSION OF TEMPERATURE CORRECTIONS (DEFAULT)
c*** when using the temperature correction procedurte:
c    this option should be switched off within steal
      NOTEMP=.TRUE.
C*** NO TEMPERATURE PLOT (DEFAILT)
      TPLOT=.FALSE.
     
!      OPEN (1,FILE='CARDS')
      OPEN (1,FILE='cards.inp')
      REWIND 1
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)

      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
               PRINT *,' STEAL: LINE BLANKETING OPTION OFF'
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
               PRINT *,' STEAL: LINE BLANKETING OPTION ON'
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
               PRINT *,' STEAL: LINE BLANKETING OPTION ON (print on)'
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
               PRINT *,' STEAL: LB TABLE INPUT'
            ELSE
               PRINT *,' NO VALID OPTION'
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:6) .EQ. 'NEWWRC' ) THEN
C                         ======
            DECODE (80,9,KARTE) XL
    9       FORMAT (7X,F10.0)
            NEWWRC=NINT(XL)
            GOTO 8
            ENDIF
      IF (KARTE(:9) .EQ. 'PLOT FLUX') THEN
C                         =========
            IPLOTF=1
            GOTO 8
            ENDIF
      IF ( KARTE(:6) .EQ. 'PROLIB' ) THEN
C                          ======
         DECODE(80,13,KARTE) TEFFE,GRAD,ALDMDT,VINF,BET,Y0
   13    FORMAT(10X,6F10.3)
            IPLOTF=1
            PROLIB=.TRUE.
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT RATE' ) THEN
                          !===========
            DECODE (80,1,KARTE) XL,XFRRA,XTORA

!0        1        2         3         4         5
!23456789 12345678 123456789 123456789 123456789 1234567890
!RINT RATE IIIIIIIII         IIIIIIIIII        IIIIIIIIII
    1       FORMAT (11X,F10.0,10X,F10.0,8X,F10.0)
            LSRAT=NINT(XL)
            IFRRA=NINT(XFRRA)
            ITORA=NINT(XTORA)
            IF (LSRAT.EQ.0) LSRAT=1
            IF (IFRRA.GT.ITORA) ITORA=IFRRA
            PRINT '("DECSTE: RATE =",I5,":",I5,":",I5)', IFRRA,LSRAT,ITORA
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT HIST' ) THEN
C                          ==========
            IHIST=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT FLUX' ) THEN
C                          ==========
            IFLUX=1
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT DATOM' ) THEN
C                          ===========
            IDAT=1
            GOTO 8
            ENDIF
      IF (KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                         =========
            DECODE (80,3,KARTE) XL
    3       FORMAT (9X,F10.0)
            LSPOP=NINT(XL)
            IF (LSPOP.EQ.0) LSPOP=1
            PRINT *,' PRINT POP ',LSPOP,' DECODED'
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT EXPO' ) THEN
C                          ==========
            DECODE (80,1,KARTE) XL
            LSEXPO=NINT(XL)
            IF (LSEXPO .EQ. 0) LSEXPO=1
            PRINT *,' PRINT EXPO ',LSEXPO,' DECODED'
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'JOBMAX=' ) THEN
C                         =======
            DECODE (80,11,KARTE) XL
   11       FORMAT (7X,F10.9)
            JOBMAX=NINT(XL)
!      IF (JOBMAX .GE. 1000) THEN
!            PRINT *,'JOBMAX .GT. 1000: INVALID OPTION'
!            STOP 'ERROR'
!            ENDIF
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'EPSILON' ) THEN
C                         =======
            DECODE (80,14,KARTE) EPSILON
   14       FORMAT (8X,F10.0)
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'REDUCE' ) THEN
C                         ======
            DECODE (80,17,KARTE) REDUCE
   17       FORMAT (7X,F10.0)
            IF (REDUCE.EQ..0) REDUCE=.5
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'PLOT POP') THEN
C                         ========
            IF (NPLOT+1 .EQ. N) GOTO 8
            DECODE (30,5,KARTE) L1,L2,L3
    5 FORMAT (9X,3(I2,1X))
            IF (L1 .LE. 0) GOTO 8
            NPLOT=NPLOT+1
            ENCODE (8,2,LEVELPL(NPLOT)) L1,L2,L3
    2       FORMAT (3I2)
            LEVELPL(NPLOT+1)=0
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT CCORE') THEN
C                          ===========
            IPRICC=1
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT LCORE') THEN
C                          ===========
            DECODE (80,7,KARTE) XLC
    7       FORMAT (11X,F10.0)
            IPRILC=NINT(XLC)
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'GAMMA=') THEN
C                         ======
            MGC=MGC+1
            MGR=MGR+1
            MGL=MGL+1
            IF (MGC .GT. 10 .OR. MGR .GT. 10 .OR. MGL .GT. 10) THEN
               write (6,*) ('TOO MANY GAMMA OPTIONS')
               STOP 'ERROR'
               ENDIF
            DECODE (80,12,KARTE) GAMMA, ANG
   12       FORMAT (6X,F11.0,9X,F10.0)
            NGAMR(MGR)=ANG
            NGAML(MGL)=ANG
            AGAMR(MGR)=GAMMA
            AGAML(MGL)=GAMMA
            GOTO 8
            ENDIF







    6       FORMAT (7X,F10.0,9X,F10.0)
      IF (KARTE(:6) .EQ. 'GAMMAR' ) THEN
C                         ======
            MGR=MGR+1
            IF (MGR .GT. 10) STOP 'ERROR'
            DECODE (80,6,KARTE) GAMMAR,ANGR
            NGAMR(MGR)=ANGR
            AGAMR(MGR)=GAMMAR
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'GAMMAL' ) THEN
C                         ======
            MGL=MGL+1
            IF (MGL .GT. 10) STOP 'ERROR'
            DECODE (80,6,KARTE) GAMMAL,ANGL
            NGAML(MGL)=ANGL
            AGAML(MGL)=GAMMAL
            GOTO 8
            ENDIF
      IF ( KARTE(:4) .EQ. 'LINE' ) THEN
C                          ====
            DECODE (80,41,KARTE) IND1
   41       FORMAT(4X,I3)

            IF (KARTE(9:10) .EQ. 'TO') THEN

!               CALL SYSTEM('grep LINE'//' '//datom_nlte//' '//'> temp.out')

!               IND2 = NUM_OF_LINES('temp.out')
               IND2 = lastind_nlte

!               CALL SYSTEM('rm temp.out')

            ELSE

               IND2 = IND1

            ENDIF

            IF (IND2 .GE. IND1) THEN
                  INC=1
                  ELSE
                  INC=-1
                  ENDIF
            DO 45 IND=IND1,IND2,INC
            IF (IND .GT. LASTIND) THEN
               PRINT *,' ERROR STOP IN DECSTE:  IND  .GT. LASTIND'
               write (6,*) (' IND  .GT. LASTIND')
               STOP 'ERROR'
               ENDIF
            LINE(IND)=.TRUE.
   45       CONTINUE
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'OLDSTART') THEN
C                         ========
            OLDSTART=.TRUE.
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                         ========
         NOTEMP=.TRUE.
         GOTO 8
         ENDIF
      IF (KARTE(:6).EQ.'PLOT T') THEN
C                       ======
            DECODE (80,15,KARTE) YMAX
   15       FORMAT (22X,F6.0)
            TPLOT=.TRUE.
         GOTO 8
            ENDIF
     
      GOTO 8

   66 CONTINUE

      CLOSE(1)

      return

      end subroutine

      subroutine decon(LSOPA,LSINT,IFLUX,JOBMAX,LPRIH,LPHNU,LPRIV,TEFF,LBLANK)

!     DECODING INPUT OPTIONS, CALLED FROM WRCONT

      implicit real*8(a - h, o - z)
      
      LOGICAL NOTEMP
      CHARACTER KARTE*80
     
C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      LPRIH=-1
      LPHNU=-1
      LPRIV=-1
      IFLUX=-1
      JOBMAX=-1
      LBLANK=0
      TEFF=-1.0d0
      NOTEMP=.FALSE.
     
!      OPEN (1,FILE='CARDS')
      OPEN (1,FILE='cards.inp')
      REWIND 1
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
               PRINT *,' WRCONT: LINE BLANKETING OPTION OFF'
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
               PRINT *,' WRCONT: LINE BLANKETING OPTION ON'
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
               PRINT *,' WRCONT: LINE BLANKETING OPTION ON (print on)'
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
               PRINT *,' WRCONT: LB TABLE INPUT'
            ELSE
               PRINT *,' NO VALID OPTIN DECODED'
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:10) .EQ. 'PRINT FLUX') THEN
C                          ==========
            IFLUX=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT INT ') THEN
C                          ==========
            read (karte,7) XL
    7       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN
C                          ==========
            read (karte,7) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT F(R)') THEN
C                          ==========
            read (karte,7) XL
            LPRIH=NINT(XL)
            IF (LPRIH.EQ.0) LPRIH=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT HNU ') THEN
C                          ==========
            read (karte,7) XL
            LPHNU=NINT(XL)
            IF (LPHNU.EQ.0) LPHNU=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT VLJ ') THEN
C                          ==========
            read (karte,7) XL
            LPRIV=NINT(XL)
            IF (LPRIV.EQ.0) LPRIV=1
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'JOBMAX=') THEN
C                         =======
            read (karte,7) XL
    3       FORMAT (7X,F10.9)
            JOBMAX=NINT(XL)
!      IF (JOBMAX .GE. 1000) THEN
!            PRINT *,'JOBMAX .GT. 1000: INVALID OPTION'
!            STOP 'ERROR'
!            ENDIF
            GOTO 8
            ENDIF
      IF (KARTE(:5) .EQ. 'TEFF=') THEN
C                         =====
            read (karte,5) TEFF
    5       FORMAT (5X,F10.0)
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                         ========
            NOTEMP=.TRUE.
            GOTO 8
            ENDIF
     
      GOTO 8

   66 CONTINUE

      close(1)

      return
     
      end subroutine

      subroutine decomo(LSOPA, LSINT, LBLANK)

!     DECODING INPUT OPTIONS FOR PROGRAM "COMO"

      implicit real*8(a-h,o-z)
     
      CHARACTER KARTE*80
!      OPEN (1,FILE='CARDS')
      OPEN (1,FILE='cards.inp')
      REWIND 1
      LSOPA=-1
      LSINT=-1
      LBLANK=0
     
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN

            PRINT *,' LINE BLANKETING OPTION DECODED COMO'
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:10) .EQ. 'PRINT INT ') THEN

            DECODE (80,7,KARTE) XL
    7       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN

            DECODE (80,7,KARTE) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            ENDIF
      GOTO 8
   66 CONTINUE

      CLOSE(1)

      return

      end subroutine

      subroutine decetl(LSOPA, LSINT, VDOP, LINE, NLINE, LINEKEY, LASTIND, LBLANK)

      use vardatom_nlte
      use file_operations

!     DECODES INPUT CARDS FOR MAIN PROGRAM "ETL"

      implicit none
      integer,intent(in   ) :: LASTIND
      real*8, intent(  out) :: VDOP
      integer,intent(  out) :: LSOPA,LSINT,NLINE,LBLANK
      logical,intent(inout) :: LINEKEY(LASTIND)
      character*7,intent(inout) :: LINE(LASTIND)
      CHARACTER :: KARTE*80
      LOGICAL :: ETLKEY
      real*8  :: XL
      integer :: IND,IND1,IND2,INC
      
!      OPEN (1,FILE='CARDS')
      OPEN (1,FILE='cards.inp')
      LSOPA=-1
      LSINT=-1
      NLINE=0
      ETLKEY=.FALSE.
      VDOP=.0
      LBLANK=0
     
    6 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            PRINT *,' ETL: LB OPTION DECODED  -  not active in this rout
     $ine'
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
            ENDIF
            GOTO 6
            ENDIF

      IF ( KARTE(:10) .EQ. 'PRINT INTL') THEN
C                           ==========
            DECODE (80,8,KARTE) XL
    8       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            GOTO 6
            ENDIF

      IF ( KARTE(:10) .EQ. 'PRINT OPAL') THEN
C                           ==========
            DECODE (80,8,KARTE) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            GOTO 6
      ELSEIF (KARTE(:5) .EQ. 'VDOP=' ) THEN
C                         =====
            DECODE (80,11,KARTE) VDOP
   11       FORMAT (5X,F10.0)
            GOTO 6
      ELSEIF ( KARTE(:4) .EQ. 'LINE' ) THEN

        read (KARTE,'(4X,I3)') IND1

        IF (KARTE(9:10) .EQ. 'TO') THEN

!           CALL SYSTEM('grep LINE'//' '//datom_nlte//' '//'> temp.out')

!           IND2 = NUM_OF_LINES('temp.out')
           IND2 = lastind_nlte

!           CALL SYSTEM('rm temp.out')

        ELSE

           IND2 = IND1

        ENDIF

        IF (IND2 .GE. IND1) THEN

          INC=1

        ELSE

          INC=-1

        ENDIF

        PRINT '("decetl:LINE ",i0," TO ",i0, ", inc=",i0)', IND1, IND2, INC

        DO IND = IND1, IND2, INC

          NLINE = NLINE + 1

          IF (NLINE .GT. LASTIND) THEN
            PRINT *,' ERROR STOP IN DECETL: NLINE .GT. LASTIND'
            ! CALL REMARK ('NLINE .GT. LASTIND')
            STOP 'ERROR'
          ENDIF
          write (LINE(NLINE),'(4HLINE,I3)') ind
          LINEKEY(NLINE)=ETLKEY
        ENDDO
        GOTO 6
      ENDIF
      IF ( KARTE(:6) .EQ. 'FORMAL' ) THEN
C                          ======
            ETLKEY=.FALSE.
            GOTO 6
            ENDIF
      IF ( KARTE(:3) .EQ. 'ETL' ) THEN
C                          ===
            ETLKEY=.TRUE.
            GOTO 6
            ENDIF
      GOTO 6
66    CONTINUE

      close(1)

      return

      end subroutine

      end module
