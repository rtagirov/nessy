      module MOD_DATOM_M

      contains

      subroutine DATOM_M(NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                   EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $                   INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                   ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $                   NLAST,WAVARR,SIGARR,NFDIM)
C2     DATOM_M
C3     | RDCSARR

C***  Changes by Margit Haberreiter, May 20, 2002
C*** CALLED BY COMO, ETL, STEAL, WRCONT, WRSTART, FIOSS8
C*******************************************************************************
CMH   DATOM_M: CODE NUMBER OF ELEMETS CHANGED
C***  READS ATOMIC DATA FROM TAPE4=DATOM  **************************************
C***  CODE NUMBER OF ELEMENT HELIUM    (HE): COMPONENT 1 OF VECTOR 'KODAT'
C***                         HYDROGEN  (H )            2
C**********************************************************
C***                         LITHIUM   (Li )           3
C***                         BERYLLIUM (Be )           4
C***                         BOR       (B )            5
C*********************************************************
C***                         CARBON    (C )            6
C**********************************************************
C***                         NITROGEN  (N )            7
C***                         OXYGEN    (O )            8
C***                         FLOUR     (FL )           9
C***                         NEON      (NE )          10
C**********************************************************
C***                         SODIUM    (Na)           11
C***                         MAGNESIUM (Mg)           12 
C***                         ALUMINIUM (Al)           13
C***                         SILICON   (Si)           14
C**********************************************************
C***                         PHOSPHOR  (P )           15
C**********************************************************
C***                         SULFUR    (S )           16
C**********************************************************
C***                         CHLOR     (CL)           17
C***                         ARGON     (AR)           18
C**********************************************************
C***                         POTASSIUM (K )           19
C***                         CALCIUM   (Ca)           20
C**********************************************************
C***                         SCANDIUM  (SC)           21
C***                         TITAN     (TI)           22
C***                         VANADIUM  (V )           23
C***                         CHROM     (CR)           24
C***                         MANGAN    (MN)           25
C**********************************************************
C***                         IRON      (Fe)           26
C**********************************************************
C***                         COBALT    (CO)           27
C**********************************************************
C***                         NICKEL    (Ni)           28
C**********************************************************
C***                         COPPER    (CU)           29
C***                         ZINK      (ZN)           30
C**********************************************************
C*** changed by Margit Haberreiter, May 2002
!    changed by Rinat Tagirov, November 2016
!    Getting rid of NDIM, MAXATOM and MAXIND
C*******************************************************************************
      use MOD_RDCSARR
      use MOD_ERROR

      implicit none

      !public variables - out

      integer, intent(out) :: N, NATOM

      integer, intent(out), allocatable, dimension(:) :: KODAT, NOM

      integer, intent(out), allocatable, dimension(:) :: NFIRST, NLAST

      integer,intent(out):: INDLOW, INDNUP
      integer,intent(out):: LASTIND,NCHARG,MAINQN
      real*8,intent(out) :: ALTESUM,ATMASS,ALPHA,COCO
      character*8,intent(out) :: agaunt
      real*8,intent(out) ::  EINST,EION, ELEVEL,SEXPO,WEIGHT,STAGE
      real*8,intent(out) ::  SIGARR, WAVARR ! get changed in RDCSARR
      !public variables - in
      integer,intent(in) :: MAXATOM,MAXIND
      integer,intent(in) :: NDIM,NFDIM
      !private variables
      integer I,IND,IRANGE,IECHO,J,LEV,LEVSEQ,LOW,LOWP
      integer MQN,NA,NCHG,NUP,NW
      real*8  ALPLOW,AUPLOW
      real*8  ASUM,CO1,CO2,CO3,CO4,COEFF1,COEFF2
      character*8 :: AGLOW
      real*8  E,ELEV,F,GFG,SLOW,SIGMA
      !constants
      real*8,parameter ::  ONE = 1.D+0
      !strings
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(NDIM),LEVUP,LEVLOW,Lread
      CHARACTER*10 ELEMENT(MAXATOM),LINECA
      CHARACTER*4 CEY,KEYCOL(NDIM,NDIM)
      CHARACTER*3 KRUDI
      CHARACTER*2 SYMBOL(MAXATOM) 
      !dimensions
!      DIMENSION NCHARG(NDIM), WEIGHT(NDIM), ELEVEL(NDIM)
!      DIMENSION EION(NDIM),MAINQN(NDIM),EINST(NDIM,NDIM)

!      DIMENSION ALPHA(NDIM),SEXPO(NDIM),AGAUNT(NDIM)
!      DIMENSION NOM(NDIM)
!      DIMENSION COCO(NDIM,NDIM,4)
!      DIMENSION ALTESUM(4,NDIM)
!      DIMENSION KODAT(MAXATOM),ATMASS(MAXATOM),STAGE(MAXATOM)
!      DIMENSION NFIRST(MAXATOM),NLAST(MAXATOM)
!      DIMENSION INDNUP(MAXIND),INDLOW(MAXIND)
!      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)

      CHARACTER*(*),parameter ::
     $FORMAT_LEVEL='(12X,A10,1X,I2,1X,I4,2F10.0,1X,I2)',
     $FORMAT_ELEMENT='(12X,A10,2X,A2,4X,F6.2,3X,F5.0)',
     $FORMAT_LTESUM='(10X,A10,1X,A8,1X,G9.0,1X,F7.0,1X,F7.0)'
      ! ------------------- Begin init -------------------------
      KODAT(:) = 0
      NATOM=0
      N=0
      LEVSEQ=0
      ALTESUM(1,:)   =-one
      COCO(:,:,:)    =0d0
      KEYCOL(:,:)    ='    '
      EINST (:,:)    =-one
      OPEN (4,FILE='DATOM',status='old',readonly)
      IECHO=0
      ! ------------- end init ------------------------------
    1 READ(4,'(A)',END=66) KARTE
c      IF (IECHO.EQ.1) PRINT *,KARTE
!    2 FORMAT(A)
      IF (KARTE(:1) .EQ. '*' ) GOTO 1     ! comment
      IF (KARTE(:4)  .EQ. 'ECHO') GOTO 3    ! outdated
      IF (KARTE(:10) .EQ. 'ELEMENT   ') GOTO 5
      IF (KARTE(:10) .EQ. 'LEVEL     ' ) GOTO 10
      IF (KARTE(:10) .EQ. 'LINE      ' ) GOTO 20
      IF (KARTE(:10) .EQ. 'CONTINUUM ' ) GOTO 30
      IF (KARTE(:10) .EQ. 'LTESUM    ') GOTO 40
c     ignore dielectronic option
      IF (KARTE(:10) .EQ. 'DIELECREC ') GOTO 1  ! ignored
      IF (KARTE(:10) .EQ. 'DRTRANSIT ') GOTO 1  ! ignored
      print *, 'UNRECOGNIZED DATA INPUT IN DATOM_M: '
      print *,'KARTE="'//KARTE//'"'
      CALL ERROR('UNRECOGNIZED DATA INPUT IN DATOM_M: "'//KARTE//'"')
    3 IECHO=1
      GOTO 1
C***  ELEMENTS ---------------------------------------------------------
    5 NATOM=NATOM+1
      if (natom .gt. 30) stop 'datom_m natom > 30'
      IF (NATOM .GT. MAXATOM) THEN
          print *, 'DATOM: MORE ELEMENTS THAN ALLOWED'
          STOP 'ERROR'
          ENDIF
      LEVSEQ=0
      read (KARTE,FORMAT_ELEMENT) ELEMENT(NATOM),SYMBOL(NATOM),
     $                    ATMASS(NATOM),STAGE(NATOM)
!    9 FORMAT (12X,A10,2X,A2,4X,F6.2,3X,F5.0)
C***  MODEL ATOM OF 'HELIUM' DECODED
      IF ( (ELEMENT(NATOM) .EQ. 'HELIUM    ') .AND.
     $     ( SYMBOL(NATOM) .EQ. 'HE' .or. SYMBOL(NATOM) .EQ. 'He' )
     $    ) THEN   
      KODAT(1)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'HYDROGEN  ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'H ')) THEN
                                     KODAT(2)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'LITHIU    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Li')) THEN
                                     KODAT(3)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'BERRYL    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Be')) THEN
                                     KODAT(4)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'BOR       ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'B ')) THEN
                                     KODAT(5)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'CARBON    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'C ')) THEN
                                     KODAT(6)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'NITROG  ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'N ')) THEN
                                     KODAT(7)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'OXYGEN    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'O ')) THEN
                                     KODAT(8)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'FLOUR    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'F ')) THEN
                                     KODAT(9)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'NEON     ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Ne')) THEN
                                     KODAT(10)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'NATRIUM   ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Na')) THEN
                                     KODAT(11)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'MAGNES    ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Mg').or.(SYMBOL(NATOM) .EQ. 'MG'))) THEN
                                     KODAT(12)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'ALUMIN    ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Al').or.(SYMBOL(NATOM) .EQ. 'AL'))) THEN
                                     KODAT(13)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'SILICON   ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Si').or.(SYMBOL(NATOM) .EQ. 'SI'))) THEN
                                     KODAT(14)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'PHOSPH    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'P ')) THEN
                                    KODAT(15)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'SULPHUR   ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'S ')) THEN
                                     KODAT(16)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'CHLOR     ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Cl')) THEN
                                    KODAT(17)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'ARGON    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Ar')) THEN
                                    KODAT(18)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'POTASS    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'K ')) THEN
                                     KODAT(19)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'CALCIUM   ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Ca').or.(SYMBOL(NATOM) .EQ. 'CA'))) THEN
                                     KODAT(20)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'SCANDI    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Sc')) THEN
                                    KODAT(21)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'TITAN    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Ti')) THEN
                                    KODAT(22)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'VANADI    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'V ')) THEN
                                    KODAT(23)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'CHROM    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Cr')) THEN
                                    KODAT(24)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'MANGAN    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Mn')) THEN
                                    KODAT(25)=NATOM
CMH  MODEL ATOM OF "IRON" DECODED
      ELSEIF ((ELEMENT(NATOM) .EQ. 'IRON      ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Fe').or.(SYMBOL(NATOM) .EQ. 'FE'))) THEN
                                     KODAT(26)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'COBALT    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Co')) THEN
                                    KODAT(27)=NATOM
CMH  MODEL ATOM OF "NICKEL" DECODED
      ELSEIF ((ELEMENT(NATOM) .EQ. 'NICKEL    ') .AND.
     $   ((SYMBOL(NATOM) .EQ. 'Ni').or.(SYMBOL(NATOM) .EQ. 'NI'))) THEN
                                     KODAT(28)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'COPPER    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Cu')) THEN
                                    KODAT(29)=NATOM
      ELSEIF ((ELEMENT(NATOM) .EQ. 'ZINC    ') .AND.
     $    (SYMBOL(NATOM) .EQ. 'Zn')) THEN
                                    KODAT(30)=NATOM
      ELSE
          print *, 'UNKNOWN ELEMENT DECODED'
          print *,KARTE
          STOP 'ERROR'
      ENDIF
      GOTO 1
     
C***  LEVELS -----------------------------------------------------------
   10 N=N+1
      if (n.gt.ndim) then
         print *,N,NDIM
         pause' N gt ndim'
         stop
         endif
      IF (LEVSEQ.NE.0) THEN
          print *, 'DATOM: LEVEL CARD OUT OF SEQUENCE'
          STOP 'ERROR'
          ENDIF
      IF(N.GT.NDIM) THEN
          print *, 'DATOM : DIMENSION OVERFLOW'
          STOP 'ERROR'
          ENDIF
      IF (NATOM .NE. 0) NOM(N)=NATOM
c      DECODE(80,11,KARTE)   LEVEL(N),NCHARG(N),NW,ELEVEL(N),E,MAINQN(N)

      nchg = 0
      nw= 0
      elev=0.
      e=0.
      mqn=0
      read(KARTE,FORMAT_LEVEL) lread,nchg,nw,elev,e,mqn
      LEVEL(N)=lread
      NCHARG(N)=nchg
      ELEVEL(N)=elev
      MAINQN(N)=mqn
CMH   IF HMINUS READ THEN SET IELHM = 1
c     PRINT *,'DATOM_M: SET IELHM =1 '
C     IF (LREAD .EQ. 'H MINUS..1') THEN 
C       IELHM=1.
C     PRINT *, 'DATOM_M ',IELHM
C     PAUSE
C     ENDIF
      WEIGHT(N)=FLOAT(NW)
      IF (ELEVEL(N).EQ..0 .AND. MAINQN(N).LE.1) EION(N)=E
C      PRINT *,'1. DATOM_M: ',EION(N)!, 1.8/EION(N)
      IF (ELEVEL(N).NE..0 .AND. MAINQN(N).LE.1 .AND. (E.NE.0.))
     $   EION(N)=E-ELEVEL(N) 
CMHpr      PRINT *,'2. DATOM_M: ',E,EION(N)!, 1.8/EION(N)
      GOTO 1
     
C***  LINE TRANSITIONS  ------------------------------------------------
   20 CONTINUE
      READ (KARTE,21) LEVUP,LEVLOW,LINECA,KRUDI,CEY,CO1,CO2,CO3,CO4
   21 FORMAT(10X,A10,2X,2A10,2X,A3,A4,1X,4G7.0)
      IF (LINECA(1:1).EQ.'F') THEN
         READ (LINECA(3:10),'(F8.4)') AUPLOW
         ELSE
         READ (LINECA(1:10),'(E10.4)') AUPLOW
         ENDIF
      LEVSEQ=1
C***  FIND UPPER INDEX
      DO 22 J=1,N
      NUP=J
      IF (LEVEL(J).EQ.LEVUP ) GOTO 23
   22 CONTINUE
      print *,levup
      print *, 'UPPER LINE LEVEL NOT FOUND'
      STOP 'ERROR'
C***  FIND LOWER INDEX
   23 CONTINUE
      DO 24 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 25
   24 CONTINUE
      print *,levlow
      print *,  'LOWER LINE LEVEL NOT FOUND'
      STOP 'ERROR'
   25 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            print *,  'LINE BETWEEN DIFFERENT ELEMENTS'
            STOP 'ERROR'
            ENDIF
         ENDIF
      IF (NCHARG(NUP) .NE. NCHARG(LOW)) THEN
            print *,  'LINE BETWEEN DIFFERENT IONIZATION STAGES'
            STOP 'ERROR'
            ENDIF
      IF (NUP.LE.LOW) THEN
            print *,  'LINE TRANSITION INDICES WRONG'
            STOP 'ERROR'
            ENDIF
      IF (LINECA(1:1).EQ.'F') THEN
         IF (LINECA(1:2).EQ.'FE') GFG=AUPLOW
         IF (LINECA(1:2).EQ.'FA') GFG=WEIGHT(LOW)*AUPLOW/WEIGHT(NUP)
         AUPLOW=GFG*0.66702*(ELEVEL(NUP)-ELEVEL(LOW))**2
         ENDIF
c Potsdam convention:
      IF (AUPLOW.lt.0.) then
         GFG=-AUPLOW
         AUPLOW=GFG*0.66702*(ELEVEL(NUP)-ELEVEL(LOW))**2
c         print *,' Potsdam-f: ',low,nup,Auplow
         ENDIF
      EINST(NUP,LOW) = AUPLOW
      KEYCOL(NUP,LOW)=CEY
      IF ((CEY(:3).EQ.'BFK') .OR. (CEY(:2).EQ.'BK')) THEN
            COCO(NUP,LOW,1)=CO1
            COCO(NUP,LOW,2)=CO2
            COCO(NUP,LOW,3)=CO3
            COCO(NUP,LOW,4)=CO4
      ENDIF
      IF ((CEY(:3).EQ.'NIT') .OR. (CEY(:3).EQ.'CAR')) THEN
            COCO(NUP,LOW,1)=CO1
            COCO(NUP,LOW,2)=0.
            COCO(NUP,LOW,3)=0.
            COCO(NUP,LOW,4)=0.
      ENDIF
C***  RUDIMENTAL TRANSITIONS ARE MARKED BY -2. IN THE TRANSPOSED
C***  MATRIX ELEMENT  EINST(LOW,NUP)
      IF (KRUDI.NE.'   ') EINST(LOW,NUP)=-2.
      GOTO 1
     
C***  CONTINUUM TRANSITIONS    -----------------------------------------
   30 DECODE (80,31,KARTE)  LEVLOW, SIGMA,ALPLOW,SLOW,AGLOW
   31 FORMAT (10X,A10,1X,3G10.0,1X,A8,1X)
C   31 FORMAT(10X,A10,2X,2A10,2X,A8,1X,A4,1X,A4)
         
         LEVSEQ=1

C***  FIND LOWER INDEX
      DO 34 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 35
   34 CONTINUE
      print *,  'LOWER CONTINUUM LEVEL NOT FOUND'
      STOP 'ERROR'
   35 CONTINUE
C***  FIND UPPER INDEX

      LOWP=LOW+1
      DO 32   J=LOWP,N
      NUP=J
c     $,'ncharg(low)',NCHARG(low),'ncharg(up)',ncharg(nup)

      IF (NCHARG(LOW)+1 .EQ. NCHARG(NUP)) GOTO 33
   32 CONTINUE
      print *,  'UPPER CONTINUUM LEVEL NOT FOUND'
      print *, ' error'
      STOP 'ERROR'
   33 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            print *,  'CONTINUUM BETWEEN DIFFERENT ELEMENTS'
            STOP 'ERROR'
            ENDIF
         ENDIF
      !***  IN CASE OF CONTINUUM CALCULATON THE PARAMETERS SIGMA,ALPLOW,SLOW, AGLOW
      !***  CAN BE MISUSED AS KEYWORDS FOR DEFINING THE TABLES FROM WHICH TO READ 
      !***  CROSS SECTIONS
      !**   number
      EINST(LOW,NUP)=SIGMA
      !*** number
      ALPHA(LOW)=ALPLOW
      !*** number
      SEXPO(LOW)=SLOW
      !*** string der Laenge 8 (E.G.'TABLE')
      AGAUNT(LOW)=AGLOW
      GOTO 1
     
C***  SUM OF TRANSITIONS TO UPPER LEVELS WHICH ARE ASSUMED TO BE IN LTE
   40 DECODE (80,FORMAT_LTESUM,KARTE) LEVLOW,IRANGE,ASUM,COEFF1,COEFF2
!   41 FORMAT (10X,A10,1X,A8,1X,G9.0,1X,F7.0,1X,F7.0)
      LEVSEQ=1
C***  FIND LOWER INDEX
      DO 42 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW) GOTO 43
   42 CONTINUE
      print *,  'LOWER LTESUM LEVEL NOT FOUND'
      STOP 'ERROR'
   43 CONTINUE
      ALTESUM(1,LOW)=ASUM
      ALTESUM(2,LOW)=COEFF1
      ALTESUM(3,LOW)=COEFF2
      ENCODE (8,44,ALTESUM(4,LOW)) IRANGE
   44 FORMAT (A8)
      GOTO 1
     
C***  END OF INPUT DATA REACHED  ---------------------------------------
   66 CLOSE (4)
     
C***  SOME CHECKS OF THE INPUT DATA ....................
      IF(N.EQ.0) THEN
            print *,  'NO ENERGY LEVELS RECOGNIZED'
            STOP 'ERROR'
            ENDIF
C***  OLD DATOM FILE (CONTAINING ONLY A HELIUM MODEL ATOM) RECOGNIZED
      IF (NATOM .EQ. 0) THEN
         NATOM=1
         ELEMENT(1)='HYDROGEN  '
         SYMBOL(1)='H '
         ATMASS(1)=one
         STAGE(1)=2.
         KODAT(1)=1
         DO 50 I=1,N
         NOM(I)=1
   50    CONTINUE
         ENDIF
C***  ALL ELEMENTS ARE CHECKED ONE BY ONE FOR EXISTENCE OF ANY LEVEL
   53 DO 54 NA=1,NATOM
      DO 56 I=1,N
      IF (NOM(I) .EQ. NA) GOTO 54
   56 CONTINUE
      print *,  'ELEMENT WITHOUT ANY LEVEL DECODED'
      STOP 'ERROR'
   54 CONTINUE
C***  LEVELS ARE CHECKED FOR CORRECT ELEMENT MEMBERSHIP
      DO 55 J=1,N
CMH      IF (LEVEL(J)(:2) .NE. SYMBOL(NOM(J))) THEN
      IF (LEVEL(J)(:1) .NE. SYMBOL(NOM(J))(:1)) THEN
         print *,  'WRONG ELEMENT MEMBERSHIP OF LEVELS'
         STOP 'ERROR'
         ENDIF
   55 CONTINUE
C***  TRANSITIONS ARE CHECKED FOR COMPLETENESS
      DO 7 I=1,N
      DO 7 J=1,N
      IF (NOM(I) .NE. NOM(J)) GOTO 7
      IF (NCHARG(I) .NE. NCHARG(J)) GOTO 8
      IF (I.LE.J) GOTO 7
      IF (EINST(I,J) .LT. .0  ) THEN
c            print *,  'LINE TRANSITION MISSING'

c            STOP 'ERROR'
            ENDIF
      GOTO 7
    8 IF (I.GE.J) GOTO 7
C***  CHARGES MUST DIFFER BY 1
      IF (NCHARG(I)+1 .NE. NCHARG(J)) GOTO 7
C***  UPPER LEVEL MUST BE GROUND STATE OF THAT ION
      IF (NCHARG(J) .EQ. NCHARG(J-1)) GOTO 7
      IF (EINST(I,J) .LT. .0 ) THEN
            print *,  'CONTINUUM TRANSITION MISSING'
            STOP 'ERROR'
            ENDIF
    7 CONTINUE
     
C***  GENERATE VECTORS NFIRST, NLAST: FIRST AND LAST LEVEL OF EACH ELEMENT
      DO 90 NA=1,NATOM
      IF (NA .EQ. 1) THEN
          NFIRST(NA)=1
      ELSE
          NFIRST(NA)=NLAST(NA-1)+1
      ENDIF
      IF (NA .LT. NATOM) THEN
        do lev=1,N
          if (nom(lev).eq.na+1) then
            nlast(na)=lev-1
            go to 91
            endif
        enddo
        stop ' isrcheq did not work'
91      continue
c          NLAST(NA)= ISRCHEQ(N,NOM(1),1,NA+1) - 1
      ELSE
          NLAST(NA)=N
      ENDIF
   90 CONTINUE
     
C***  GENERATE VECTORS INDNUP, INDLOW: LEVEL INCICES OF THE LINES
      DO 94 IND=1,MAXIND
      INDNUP(IND)=0
      INDLOW(IND)=0
   94 CONTINUE
      IND=0
      DO 95 NUP=2,N
      DO 95 LOW=1,NUP-1
      IF ((NCHARG(LOW) .NE. NCHARG(NUP)) .OR. (NOM(LOW) .NE. NOM(NUP)))
     $   GOTO 95
      IND=IND+1
      INDNUP(IND)=NUP
      INDLOW(IND)=LOW
   95 CONTINUE
      LASTIND=IND
C***  ERROR STOP
      IF (LASTIND .GT. MAXIND) THEN
         print *, 'LASTIND .GT. MAXIND'
      STOP 'ERROR'
      ENDIF
     
C***  ASSIGNMENT OF IONIZATION ENERGIES (INPUT VALUES OF GROUND STATE)
C***  TO ALL LEVELS OF THE CORRESPONDING ELEMENT
      DO 13 I=1,N
      IF (EION(I) .NE. 0.) THEN
         DO 14 J=1,N
         IF ((NOM(J) .EQ. NOM(I)) .AND. (NCHARG(J) .EQ. NCHARG(I)))
     $      EION(J)=EION(I)
   14    CONTINUE
         ENDIF
   13 CONTINUE
     
C***  IF MAIN QUANTUM NUMBER IS GIVEN, LEVEL ENERGIES ARE COMPUTED BY
C***  RYDBERG FORMULA
      J=0
      DO 4 J=1,N
      IF (MAINQN(J).LE.1) GOTO 4
      F=FLOAT(MAINQN(J))
      ELEVEL(J)=(one-one/F/F)*EION(J)
    4 CONTINUE

C***  IF AGAUNT(LOW) EQ 'TABLE' READ CROSS SECTION FROM TABLE
      J=0 
      DO J=1,N
        ! LOW=J
        IF (AGAUNT(J) .EQ. 'TABLE') THEN 
          ! print *,j,n 
          ! read in WAVARR and SIGARR from the File LEVEL(J)
          call RDCSARR(LEVEL,J,N,WAVARR,SIGARR,NDIM,NFDIM)
          ! PRINT *,J,N
        END IF
      END DO
C********************************************************************

c      do i = 1, ndim

c         print*, 'LEVELS CHECK', i, LEVEL(i)

c      enddo

c      do i = 1, NATOM

c         print*, 'ELEMENTS CHECK', i, KODAT(i), SYMBOL(i), ELEMENT(i)

c      enddo

c      STOP

      RETURN
      END subroutine
      end module
