      MODULE MOD_DECSTAR_M

        type :: DECSTAR_PARAMS 
          integer :: POPNUM_CP
        end type
        type(DECSTAR_PARAMS) :: DECSTAR_OPT = DECSTAR_PARAMS(0)

      CONTAINS

      SUBROUTINE DECSTAR_M(MODHEAD,FM,RSTAR,VDOP,TTABLE,FAL,LBKG,XLBKG1,XLBKG2,
     $                     TPLOT,NATOM,ABXYZ,KODAT,IDAT,LBLANK,ATMEAN)
C***********************************************************************
C***  DECODES INPUT CARDS, CALLED FROM WRSTART
C***********************************************************************
      USE MOD_INITVEL
      USE FILE_OPERATIONS
      USE COMMON_BLOCK
 
      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER :: HEIGHT_DIM

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
      DIMENSION ABXYZ(NATOM),KODAT(NATOM)
      LOGICAL TTABLE, SPHERIC, TPLOT, FAL, LBKG
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
       INTEGER XLBKG1,XLBKG2
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
      DATA RSUN / 6.96D10 /

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

      OPEN (1, FILE = 'CARDS', STATUS = 'OLD'); REWIND 1

   10 READ (1,11,END=66) KARTE
   11 FORMAT (A)
      PRINT 2,KARTE
    2 FORMAT (1X,A)
     
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
               PRINT *,' NO VALID OPTIN'
               LBLANK=0
            ENDIF

            GOTO 10
            ENDIF

      IF (KARTE(:11) .EQ. 'PRINT DATOM' ) THEN
C                          ===========
        IDAT=1
        GOTO 10
      ENDIF
C***  FALP IS THE LOGICAL KEY FOR THE RISING TEMPERATURE PROFILE
      IF (KARTE(:3).EQ.'FAL') THEN
C                       =====
        FAL=.TRUE.
        GOTO 10
      ENDIF
      IF (KARTE(:4).EQ.'LBKG') THEN
C                       =====
        LBKG=.TRUE.
c       read (KARTE,*) XLBKG1,XLBKG2
        DECODE (80,100,KARTE) XLBKG1,XLBKG2
  100   FORMAT (5X,I4,1X,I5)
        PRINT *,'LBKG= ',LBKG,XLBKG1,XLBKG2
        IF ((XLBKG1 .eq. 0) .or. (XLBKG1 .eq. 0) .or.
     $              (xlbkg2 .lt. xlbkg1)) THEN
          write (6,*) 'something wrong with lbkg-wavelength range'
          write (6,*) 'decstar_m: LBKG= ',LBKG,XLBKG1,XLBKG2
          print *, 'possibly wrong format lbkg-wavelength range'
          print *, 'decstar_m: LBKG= ',LBKG,XLBKG1,XLBKG2
          pause
        ENDIF
        GOTO 10
      ENDIF
      IF (KARTE(:5).EQ.'TABLE') THEN
C                       =====
            TTABLE=.TRUE.
            GOTO 10
            ENDIF
      IF (KARTE(:5).EQ.'SPHER') THEN
C                       =====
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
      RSTAR=RSTAR*RSUN
C***  INITIALISATION OF THE VELOCITY-FIELD PARAMETERS
      IF (RPAR(20:24).NE.'TABLE') THEN

!        RINAT TAGIROV:
!        This bit is necessary to correct for the 
!        inconsistency of the RMAX parameter given in the CARDS file.
!        Here we change the RMAX value, so that it becomes consistent 
!        with the height grid from the FAL_VD file.
!        Without this correction the radial scale in the velocity law becomes 
!        inconsistent with RMAX and the velocity boundary values calculated in the code
!        no longer match those indicated in the CARDS file.
!        See INITVEL.FOR, WRVEL.FOR and GEOMESH.FOR for details.
!*************************************************************************************
         HEIGHT_DIM = NUM_OF_LINES(fal_mod_file)

         IF (ALLOCATED(HEIGHT)) deallocate(HEIGHT)

         allocate(height(height_dim))

         HEIGHT = READ_ATM_MOD(fal_mod_file, '1')

         RMAX = 1.0 + MAXVAL(HEIGHT) * 1.0D+5 / RSTAR
!*************************************************************************************

         CALL INITVEL(RMAX,TEFF,GLOG,RSTAR,XMASS)

      ENDIF

      RETURN

      END SUBROUTINE

      END MODULE
