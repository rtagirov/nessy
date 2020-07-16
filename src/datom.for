      module mod_datom

      contains

      subroutine datom(mode,
     $                 N,
     $                 LEVEL,
     $                 NCHARG,
     $                 WEIGHT,
     $                 ELEVEL,
     $                 EION,
     $                 MAINQN,
     $                 EINST,
     $                 ALPHA,
     $                 SEXPO,
     $                 AGAUNT,
     $                 COCO,
     $                 KEYCOL,
     $                 ALTESUM,
     $                 INDNUP, 
     $                 INDLOW,
     $                 LASTIND,
     $                 NATOM,
     $                 ELEMENT,
     $                 SYMBOL,
     $                 NOM,
     $                 KODAT,
     $                 ATMASS,
     $                 STAGE,
     $                 NFIRST,
     $                 NLAST,
     $                 WAVARR,
     $                 SIGARR,
     $                 eleatnum,
     $                 levatnum,
     $                 NFDIM)

C*** CALLED BY COMO, ETL, STEAL, WRCONT, WRSTART, FIOSS
C*******************************************************************************
CMH   DATOM: CODE NUMBER OF ELEMETS CHANGED
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
C***                         FLUOR     (FL )           9
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
C*********************************************************************************

      use utils
      use file_operations

      implicit none

      character (len = 4), intent(in)                        :: mode

      integer, intent(in)                                    :: NFDIM

      integer, intent(out)                                   :: N, NATOM, LASTIND

      integer, intent(out), allocatable, dimension(:)        :: KODAT, NOM, NFIRST, NLAST, INDLOW, INDNUP, NCHARG, MAINQN

!     eleatnum - atomic number of an element, levatnum - atomic number of a level
      integer, intent(out), allocatable, dimension(:)        :: eleatnum, levatnum

      real*8,  intent(out), allocatable, dimension(:)        :: ATMASS, ALPHA, EION, ELEVEL, SEXPO, WEIGHT, STAGE

      real*8,  intent(out), allocatable, dimension(:, :)     :: ALTESUM, EINST, SIGARR, WAVARR

      real*8,  intent(out), allocatable, dimension(:, :, :)  :: COCO

      character*2, intent(out), allocatable, dimension(:)    :: symbol

      character*4, intent(out), allocatable, dimension(:, :) :: keycol

      character*8, intent(out), allocatable, dimension(:)    :: agaunt

      character*10, intent(out), allocatable, dimension(:)   :: level, element

      !private variables
      integer ::     I, IND, IRANGE, IECHO, J, LEV, LEVSEQ, LOW, LOWP
      integer ::     MQN, NA, NCHG, NUP, NW
      real*8 ::      ALPLOW, AUPLOW
      real*8 ::      ASUM, CO1, CO2, CO3, CO4, COEFF1, COEFF2
      real*8 ::      E, ELEV, F, GFG, SLOW, SIGMA
      character*8 :: AGLOW

      integer ::                            elenum, levnum, linnum

      logical ::                            nlte_elem

      integer ::                            un

!      integer ::                            num_i_stages

!      integer, dimension(:), allocatable :: num_stage_lev

      !constants
      real*8, parameter :: ONE = 1.D+0

      !strings
      CHARACTER*80 KARTE
      CHARACTER*10 LEVUP, LEVLOW, LINECA, lread
      CHARACTER*4  CEY
      CHARACTER*3  KRUDI

      CHARACTER*(*), parameter ::
     $ FORMAT_LEVEL =   '(12X,A10,1X,I2,1X,I4,2F10.0,1X,I2)',
     $ FORMAT_ELEMENT = '(12X,A10,2X,A2,4X,F6.2,3X,F5.0)',
     $ FORMAT_LTESUM =  '(10X,A10,1X,A8,1X,G9.0,1X,F7.0,1X,F7.0)'

!     ---------------------- Array Allocation ------------------

      if (allocated(kodat))    deallocate(kodat)
      if (allocated(nom))      deallocate(nom)
      if (allocated(eleatnum)) deallocate(eleatnum)
      if (allocated(levatnum)) deallocate(levatnum)
      if (allocated(nfirst))   deallocate(nfirst)
      if (allocated(nlast))    deallocate(nlast)
      if (allocated(indlow))   deallocate(indlow)
      if (allocated(indnup))   deallocate(indnup)
      if (allocated(ncharg))   deallocate(ncharg)
      if (allocated(mainqn))   deallocate(mainqn)
      if (allocated(atmass))   deallocate(atmass)
      if (allocated(alpha))    deallocate(alpha)
      if (allocated(eion))     deallocate(eion)
      if (allocated(elevel))   deallocate(elevel)
      if (allocated(sexpo))    deallocate(sexpo)
      if (allocated(weight))   deallocate(weight)
      if (allocated(stage))    deallocate(stage)
      if (allocated(altesum))  deallocate(altesum)
      if (allocated(einst))    deallocate(einst)
      if (allocated(sigarr))   deallocate(sigarr)
      if (allocated(wavarr))   deallocate(wavarr)
      if (allocated(coco))     deallocate(coco)
      if (allocated(symbol))   deallocate(symbol)
      if (allocated(keycol))   deallocate(keycol)
      if (allocated(agaunt))   deallocate(agaunt)
      if (allocated(level))    deallocate(level)
      if (allocated(element))  deallocate(element)

      call datom_nums(mode, elenum, levnum, linnum)

      print*, mode, elenum, levnum, linnum

      allocate(kodat(30))
      allocate(nfirst(elenum))
      allocate(nlast(elenum))
      allocate(atmass(elenum))
      allocate(stage(elenum))
      allocate(symbol(elenum))
      allocate(element(elenum))
      allocate(eleatnum(elenum))

      allocate(nom(levnum))
      allocate(levatnum(levnum))
      allocate(ncharg(levnum))
      allocate(mainqn(levnum))
      allocate(alpha(levnum))
      allocate(eion(levnum))
      allocate(elevel(levnum))
      allocate(sexpo(levnum))
      allocate(weight(levnum))
      allocate(agaunt(levnum))
      allocate(level(levnum))

      allocate(coco(levnum, levnum, 4))
      allocate(keycol(levnum, levnum))
      allocate(einst(levnum, levnum))
      allocate(altesum(4, levnum))

      allocate(indlow(linnum))
      allocate(indnup(linnum))

      allocate(sigarr(levnum, NFDIM))
      allocate(wavarr(levnum, NFDIM))

      ! ------------------- Begin init -------------------------

      NATOM = 0
      N = 0
      LEVSEQ = 0

      KODAT(:) = 0
      ALTESUM(1, :) = -one
      COCO(:, :, :) = 0d0
      KEYCOL(:, :)  = '    '
      EINST (:, :)  = -one

      nfirst(:) = 0
      nlast(:) =  0

      indlow(:) = 0
      indnup(:) = 0

      nom(:) =      0
      eleatnum(:) = 0
      levatnum(:) = 0
      ncharg(:) =   0
      mainqn(:) =   0

      atmass(:) =  0.0D0
      alpha(:) =   0.0D0

      stage(:) =   0.0D0
      eion(:) =    0.0D0
      elevel(:) =  0.0D0

      sexpo(:) =   0.0D0
      weight(:) =  0.0D0

      agaunt(:) = '        '

      wavarr(:, :) = 0.0d0
      sigarr(:, :) = 0.0d0

      open(unit = 4, file = atomic_data_file, action = 'read')

      iecho = 0

      ! ------------- end init ------------------------------

    1 read(4, '(A)', end = 66) karte

      if (karte(:1) == '*' )               goto 1 ! ignore lines startring with '*'

      if     (mode == 'full') then

          if (karte(:10) == 'ELEMENT   ')  goto 5

!          if (karte(:10) == 'LEVEL     ' ) goto 10

!          if (karte(:10) == 'LINE      ' ) goto 20

!          if (karte(:10) == 'CONTINUUM ' ) goto 30

      elseif (mode == 'nlte') then

          if (karte(:10) == 'ELEMENT   ')  then

              if (index(karte, 'NLTE') /= 0) then

                  nlte_elem = .true.

                  goto 5

              else

                  nlte_elem = .false.

                  goto 1

              endif

          endif

          if (.not. nlte_elem) goto 1

!          if (karte(:10) == 'LEVEL     ') goto 10

!          if (karte(:10) == 'LINE      ') goto 20

!          if (karte(:10) == 'CONTINUUM ') goto 30

      endif

      if (karte(:10) == 'LEVEL     ') goto 10

      if (karte(:10) == 'LINE      ') goto 20

      if (karte(:10) == 'CONTINUUM ') goto 30

      call error('datom.for: mode = '//mode//'; unrecognized data input in '//
     $           atomic_data_file//':      '//karte)

    3 iecho = 1

      goto 1

!     ELEMENTS ---------------------------------------------------------
    5 NATOM = NATOM + 1

      if (natom .gt. 30) stop 'datom: natom > 30'

      LEVSEQ = 0
      read (KARTE,FORMAT_ELEMENT) ELEMENT(NATOM),SYMBOL(NATOM),ATMASS(NATOM),STAGE(NATOM)

      IF ((ELEMENT(NATOM) .EQ. 'HELIUM    ') .AND. (SYMBOL(NATOM) .EQ. 'HE' .or. SYMBOL(NATOM) .EQ. 'He')) THEN

             KODAT(1) = NATOM; eleatnum(NATOM) = 2

      ELSEIF ((ELEMENT(NATOM) .EQ. 'HYDROGEN  ') .AND. (SYMBOL(NATOM) .EQ. 'H ')) THEN

             KODAT(2) = NATOM; eleatnum(NATOM) = 1

      ELSEIF ((ELEMENT(NATOM) .EQ. 'LITHIU    ') .AND. (SYMBOL(NATOM) .EQ. 'Li')) THEN

             KODAT(3) = NATOM; eleatnum(NATOM) = 3

      ELSEIF ((ELEMENT(NATOM) .EQ. 'BERRYL    ') .AND. (SYMBOL(NATOM) .EQ. 'Be')) THEN

             KODAT(4) = NATOM; eleatnum(NATOM) = 4

      ELSEIF ((ELEMENT(NATOM) .EQ. 'BOR       ') .AND. (SYMBOL(NATOM) .EQ. 'B ')) THEN

             KODAT(5) = NATOM; eleatnum(NATOM) = 5

      ELSEIF ((ELEMENT(NATOM) .EQ. 'CARBON    ') .AND. (SYMBOL(NATOM) .EQ. 'C ')) THEN

             KODAT(6) = NATOM; eleatnum(NATOM) = 6

      ELSEIF ((ELEMENT(NATOM) .EQ. 'NITROG  ') .AND. (SYMBOL(NATOM) .EQ. 'N ')) THEN

             KODAT(7) = NATOM; eleatnum(NATOM) = 7

      ELSEIF ((ELEMENT(NATOM) .EQ. 'OXYGEN    ') .AND. (SYMBOL(NATOM) .EQ. 'O ')) THEN

             KODAT(8) = NATOM; eleatnum(NATOM) = 8

      ELSEIF ((ELEMENT(NATOM) .EQ. 'FLUOR    ') .AND. (SYMBOL(NATOM) .EQ. 'F ')) THEN

             KODAT(9) = NATOM; eleatnum(NATOM) = 9

      ELSEIF ((ELEMENT(NATOM) .EQ. 'NEON     ') .AND. (SYMBOL(NATOM) .EQ. 'Ne')) THEN

             KODAT(10) = NATOM; eleatnum(NATOM) = 10

      ELSEIF ((ELEMENT(NATOM) .EQ. 'NATRIUM   ') .AND. (SYMBOL(NATOM) .EQ. 'Na')) THEN

             KODAT(11) = NATOM; eleatnum(NATOM) = 11

      ELSEIF ((ELEMENT(NATOM) .EQ. 'MAGNES    ') .AND. ((SYMBOL(NATOM) .EQ. 'Mg') .or. (SYMBOL(NATOM) .EQ. 'MG'))) THEN

             KODAT(12) = NATOM; eleatnum(NATOM) = 12

      ELSEIF ((ELEMENT(NATOM) .EQ. 'ALUMIN    ') .AND. ((SYMBOL(NATOM) .EQ. 'Al').or.(SYMBOL(NATOM) .EQ. 'AL'))) THEN

             KODAT(13) = NATOM; eleatnum(NATOM) = 13

      ELSEIF ((ELEMENT(NATOM) .EQ. 'SILICON   ') .AND. ((SYMBOL(NATOM) .EQ. 'Si').or.(SYMBOL(NATOM) .EQ. 'SI'))) THEN

             KODAT(14) = NATOM; eleatnum(NATOM) = 14

      ELSEIF ((ELEMENT(NATOM) .EQ. 'PHOSPH    ') .AND. (SYMBOL(NATOM) .EQ. 'P ')) THEN

             KODAT(15) = NATOM; eleatnum(NATOM) = 15

      ELSEIF ((ELEMENT(NATOM) .EQ. 'SULPHUR   ') .AND. (SYMBOL(NATOM) .EQ. 'S ')) THEN

             KODAT(16) = NATOM; eleatnum(NATOM) = 16

      ELSEIF ((ELEMENT(NATOM) .EQ. 'CHLOR     ') .AND. (SYMBOL(NATOM) .EQ. 'Cl')) THEN

             KODAT(17) = NATOM; eleatnum(NATOM) = 17

      ELSEIF ((ELEMENT(NATOM) .EQ. 'ARGON    ') .AND. (SYMBOL(NATOM) .EQ. 'Ar')) THEN

             KODAT(18) = NATOM; eleatnum(NATOM) = 18

      ELSEIF ((ELEMENT(NATOM) .EQ. 'POTASS    ') .AND. (SYMBOL(NATOM) .EQ. 'K ')) THEN

             KODAT(19) = NATOM; eleatnum(NATOM) = 19

      ELSEIF ((ELEMENT(NATOM) .EQ. 'CALCIUM   ') .AND. ((SYMBOL(NATOM) .EQ. 'Ca') .or. (SYMBOL(NATOM) .EQ. 'CA'))) THEN

             KODAT(20) = NATOM; eleatnum(NATOM) = 20

      ELSEIF ((ELEMENT(NATOM) .EQ. 'SCANDI    ') .AND. (SYMBOL(NATOM) .EQ. 'Sc')) THEN

             KODAT(21) = NATOM; eleatnum(NATOM) = 21

      ELSEIF ((ELEMENT(NATOM) .EQ. 'TITAN    ') .AND. (SYMBOL(NATOM) .EQ. 'Ti')) THEN

             KODAT(22) = NATOM; eleatnum(NATOM) = 22

      ELSEIF ((ELEMENT(NATOM) .EQ. 'VANADI    ') .AND. (SYMBOL(NATOM) .EQ. 'V ')) THEN

             KODAT(23) = NATOM; eleatnum(NATOM) = 23

      ELSEIF ((ELEMENT(NATOM) .EQ. 'CHROM    ') .AND. (SYMBOL(NATOM) .EQ. 'Cr')) THEN

             KODAT(24) = NATOM; eleatnum(NATOM) = 24

      ELSEIF ((ELEMENT(NATOM) .EQ. 'MANGAN    ') .AND. (SYMBOL(NATOM) .EQ. 'Mn')) THEN

             KODAT(25) = NATOM; eleatnum(NATOM) = 25

CMH  MODEL ATOM OF "IRON" DECODED
      ELSEIF ((ELEMENT(NATOM) .EQ. 'IRON      ') .AND. ((SYMBOL(NATOM) .EQ. 'Fe') .or. (SYMBOL(NATOM) .EQ. 'FE'))) THEN

             KODAT(26) = NATOM; eleatnum(NATOM) = 26

      ELSEIF ((ELEMENT(NATOM) .EQ. 'COBALT    ') .AND. (SYMBOL(NATOM) .EQ. 'Co')) THEN

             KODAT(27) = NATOM; eleatnum(NATOM) = 27

CMH  MODEL ATOM OF "NICKEL" DECODED
      ELSEIF ((ELEMENT(NATOM) .EQ. 'NICKEL    ') .AND. ((SYMBOL(NATOM) .EQ. 'Ni') .or. (SYMBOL(NATOM) .EQ. 'NI'))) THEN

             KODAT(28) = NATOM; eleatnum(NATOM) = 28

      ELSEIF ((ELEMENT(NATOM) .EQ. 'COPPER    ') .AND. (SYMBOL(NATOM) .EQ. 'Cu')) THEN

             KODAT(29) = NATOM; eleatnum(NATOM) = 29

      ELSEIF ((ELEMENT(NATOM) .EQ. 'ZINC    ') .AND. (SYMBOL(NATOM) .EQ. 'Zn')) THEN

             KODAT(30) = NATOM; eleatnum(NATOM) = 30

      ELSE

          print*, 'UNKNOWN ELEMENT DECODED'
          print*, KARTE
          STOP    'ERROR'

      ENDIF

      goto 1

!     LEVELS -----------------------------------------------------------
   10 N = N + 1

      IF (LEVSEQ .NE. 0) THEN

          print*, 'DATOM: LEVEL CARD OUT OF SEQUENCE'

          STOP 'ERROR'

      ENDIF

      IF (NATOM .NE. 0) THEN

         NOM(N) = NATOM

         levatnum(N) = eleatnum(NATOM)

      ENDIF

      nchg = 0
      nw =   0
      elev = 0.
      e =    0.
      mqn =  0

      read(KARTE, FORMAT_LEVEL) lread, nchg, nw, elev, e, mqn

      LEVEL(N) = lread
      NCHARG(N) = nchg
      ELEVEL(N) = elev
      MAINQN(N) = mqn
      WEIGHT(N) = FLOAT(NW)

      IF (ELEVEL(N).EQ..0 .AND. MAINQN(N).LE.1) EION(N)=E

      IF (ELEVEL(N).NE..0 .AND. MAINQN(N).LE.1 .AND. (E.NE.0.))
     $   EION(N)=E-ELEVEL(N) 

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
   66 CLOSE(4)
     
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
      DO 94 IND=1,linnum
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

!     IF AGAUNT(LOW) == 'TABLE' READ CROSS SECTION FROM TABLE
      do j = 1, N

         if (AGAUNT(j) == 'TABLE') call rdcsarr(level, j, wavarr, sigarr, levnum, NFDIM)

      enddo

      return

      end subroutine

      subroutine datom_nums(mode, elenum, levnum, linnum)

      use utils
      use file_operations

      character (len = 4), intent(in) :: mode

      integer, intent(out)            :: elenum, levnum, linnum

      character (len = 100)           :: str

      logical                         :: element, level, line, continuum

      logical                         :: nlte_elem

      integer                         :: un, io

      if (mode /= 'nlte' .and. mode /= 'full')
     $stop 'datom.for: subroutine datom_nums: mode is not recognized. Abort.'

      elenum = 0
      levnum = 0
      linnum = 0

      un = 1743; open(unit = un, file = atomic_data_file, action = 'read')

      io = 0

      if     (mode == 'full') then

          do while (io == 0)

    1         read(un, '(A)', iostat = io) str

              if (io /= 0) exit

              if (str(:1) == '*' ) goto 1 ! ignore lines startring with '*'

              element =   str(:10) == 'ELEMENT   '
              level =     str(:10) == 'LEVEL     '
              line =      str(:10) == 'LINE      '
              continuum = str(:10) == 'CONTINUUM '

              if (element)   elenum = elenum + 1
              if (level)     levnum = levnum + 1
              if (line)      linnum = linnum + 1
              if (continuum) goto 1

              if (.not. element .and. .not. level .and. .not. line .and. .not. continuum)
     $        call error('datom.for: subroutine datom_nums: mode = '//mode//
     $                   '; unrecognized data input in '//atomic_data_file//':      '//str)

          enddo

      elseif (mode == 'nlte') then

          do while (io == 0)

    2         read(un, '(A)', iostat = io) str

              if (io /= 0) exit

              if (str(:1) == '*' ) goto 2 ! ignore lines startring with '*'

              element =   str(:10) == 'ELEMENT   '
              level =     str(:10) == 'LEVEL     '
              line =      str(:10) == 'LINE      '
              continuum = str(:10) == 'CONTINUUM '

              if (element) then

                  if (index(str, 'NLTE') /= 0) then

                      nlte_elem = .true.

                      elenum = elenum + 1

                      goto 2

                  else

                      nlte_elem = .false.

                      goto 2

                  endif

              endif

              if (.not. nlte_elem) goto 2

              if (level    ) levnum = levnum + 1
              if (line     ) linnum = linnum + 1
              if (continuum) goto 2

              if (.not. element .and. .not. level .and. .not. line .and. .not. continuum)
     $        call error('datom.for: subroutine datom_nums: mode = '//mode//
     $                   '; unrecognized data input in '//atomic_data_file//':      '//str)

          enddo

      endif

      close(un)

      end subroutine datom_nums


      subroutine rdcsarr(level, j, wavarr, sigarr, N, NFDIM)

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION WAVARR(N, NFDIM), SIGARR(N, NFDIM)

CMH   READS WAVENUMBER AND CROSS SECTIONS FOR EACH EXPLICIT LEVEL INTO AN ARRAY
CMH   N :     NUMBER OF LEVELS
CMH   NFDIM:  MAX NUMBER OF FREQUENCIES
CMH   WAVARR: WAVENUMBERS FOR EACH LEVEL
CMH   SIGARR: CROSS SECTIONS FOR EACH LEVEL

      integer, intent(in) :: J, N, NFDIM
      CHARACTER*10, intent(in) :: LEVEL(N)

      real*8, intent(out) :: WAVARR, SIGARR

      INTEGER IOSTATUS
      CHARACTER*10, FILENAME
      real*8 WLTH, SIGMA

      k = 1

      SIGARR(J, :) = 0.
      WAVARR(J, :) = 0.

      FILENAME = LEVEL(J)

      open(unit = 1, file = './crs/'//FILENAME, STATUS='OLD', IOSTAT=IOSTATUS, err=888, action='read')

      do while (IOSTATUS .eq. 0)

         read(unit = 1, fmt=*, IOSTAT=IOSTATUS), WLTH, SIGMA

         if (wlth .eq. 0.0) then

             print *,'RDCSARR: PROBLEM CROSS SECTION INPUT!'

             stop

         endif

         wavarr(J, K) = 1.0D8 / wlth
         sigarr(J, K) = sigma

         k = k + 1

      enddo

888   continue

      if (k .eq. 1) then

        print *,'RDCSARR: PROBLEM with ',FILENAME

        stop

      endif

      close(unit = 1)

    1 continue

      return

      end subroutine

      end module
