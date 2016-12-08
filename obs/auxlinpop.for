      MODULE LINPOP_AUX

      CONTAINS

      SUBROUTINE PRINT_NLTE_LEV(Level,
     $                          LevPopLTE,
     $                          LevPop,
     $                          DepartZwaan)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           Level

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LevPop, LevPopLTE, DepartZwaan

      CHARACTER(:), ALLOCATABLE ::          FILE_NAME, LevName

      REAL*8 ::                             RAND

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      IF (Level(1:3) .NE. 'HEI') LevName = RM_CHAR(RM_CHAR(RM_CHAR(Level, ' '), '.'), '-')

      IF (Level(1:3) .EQ. 'HEI') LevName = RM_CHAR(Level, ' ')

      IF (LevName .EQ. 'HMINUS1')   LevName = 'HMINUS'
      IF (LevName .EQ. 'ELECTRONS') LevName = 'ELECTR'

      FILE_NAME = TRIM(ADJUSTL(NLTE_DIR_1//LevName))

      CALL RANDOM_NUMBER(RAND)

      FILE_UNIT = FLOOR(100D0 + RAND * 1000D0)

      IF (.NOT. NLTE_REWRITE) THEN

         CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

         IF (LAMBDA_ITER .EQ. 0) WRITE(FILE_UNIT, '(3x,A,3x,A,8x,A,12x,A,14x,A,/)') 'LI',
     $                                                                              'DI',
     $                                                                              'LPLTE',
     $                                                                              'LP',
     $                                                                              'DZ'

      ELSE

         CALL RM_FILE(FILE_NAME, '-f')

         CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)
         
         WRITE(FILE_UNIT, '(3x,A,3x,A,8x,A,12x,A,14x,A,/)') 'LI', 'DI', 'LPLTE', 'LP', 'DZ'

      ENDIF

      DO DI = 1, DPN

         WRITE(FILE_UNIT, '(I5,2x,I3,2x,E15.7,1x,E15.7,1x,E15.7)')
     $         LAMBDA_ITER, DI, LevPopLTE(DI), LevPop(DI), DepartZwaan(DI)

      ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_NLTE_LEV


      SUBROUTINE PRINT_HYD_NLTE_TRA(LowLevel, UpLevel,
     $                              LineInd,
     $                              WAV,
     $                              VDOP,
     $                              LO,
     $                              JOLD,
     $                              SLOLD,
     $                              JNEW,
     $                              SLNEW,
     $                              PopLow, PopUp,
     $                              RadUpLow, RadLowUp,
     $                              ColUpLow, ColLowUp,
     $                              R_BRACKET)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           UpLevel, LowLevel

      INTEGER, INTENT(IN) ::                LineInd

      REAL*8, INTENT(IN) ::                 WAV, VDOP

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LO, SLNEW, JOLD, SLOLD, JNEW, R_BRACKET

      REAL*8, DIMENSION(DPN), INTENT(IN) :: PopLow, PopUp

      REAL*8, DIMENSION(DPN), INTENT(IN) :: RadUpLow, RadLowUp, ColUpLow, ColLowUp

      CHARACTER(:), ALLOCATABLE ::          FILE_NAME, TRAN_NAME, LowLev, UpLev

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      LowLev = RM_CHAR(RM_CHAR(LowLevel, ' '), '.')

      IF (LowLev .EQ. 'HMINUS1') LowLev = 'HMINUS'

      UpLev = RM_CHAR(RM_CHAR(UpLevel, ' '), '.')

      TRAN_NAME = LowLev//'_'//UpLev

      FILE_NAME = TRIM(ADJUSTL(NLTE_DIR_2//TRAN_NAME))

      FILE_UNIT = LineInd * 165

      IF (.NOT. NLTE_REWRITE) THEN

         CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

         IF (LAMBDA_ITER .EQ. 0)
     $   WRITE(FILE_UNIT, '(3x,A,7x,A,2x,A,6x,A,7x,A,7x,A,2x,A,13x,A,19x,A,20x,A,17x,A,12x,A,3(14x,A),4(13x,A),/)') 'LL',
     $                                                                                                              'UL',
     $                                                                                                              'IND',
     $                                                                                                              'WAV',
     $                                                                                                              'VDOP',
     $                                                                                                              'LI',
     $                                                                                                              'DI',
     $                                                                                                              'LOE',
     $                                                                                                              'JOL',
     $                                                                                                              'SLO',
     $                                                                                                              'JNE',
     $                                                                                                              'SLN',
     $                                                                                                              'PU',
     $                                                                                                              'PL',
     $                                                                                                              'RUL',
     $                                                                                                              'RLU',
     $                                                                                                              'CUL',
     $                                                                                                              'CLU',
     $                                                                                                              'RBR'

      ELSE

         CALL RM_FILE(FILE_NAME, '-f')

         CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

         WRITE(FILE_UNIT, '(3x,A,7x,A,2x,A,6x,A,7x,A,7x,A,2x,A,13x,A,19x,A,20x,A,17x,A,12x,A,3(14x,A),4(13x,A),/)') 'LL',
     $                                                                                                              'UL',
     $                                                                                                              'IND',
     $                                                                                                              'WAV',
     $                                                                                                              'VDOP',
     $                                                                                                              'LI',
     $                                                                                                              'DI',
     $                                                                                                              'LOE',
     $                                                                                                              'JOL',
     $                                                                                                              'SLO',
     $                                                                                                              'JNE',
     $                                                                                                              'SLN',
     $                                                                                                              'PU',
     $                                                                                                              'PL',
     $                                                                                                              'RUL',
     $                                                                                                              'RLU',
     $                                                                                                              'CUL',
     $                                                                                                              'CLU',
     $                                                                                                              'RBR'

      ENDIF

      DO DI = 1, DPN

         WRITE(FILE_UNIT, '(2(A6,2x),I3,3x,ES9.3,3x,ES7.1,2x,I5,1x,I3,1x,3(E23.15),9(1x,E15.7))')
     $         LowLev, UpLev, LineInd, WAV, VDOP, LAMBDA_ITER, DI,
     $         LO(DI),
     $         JOLD(DI),
     $         SLOLD(DI),
     $         JNEW(DI),
     $         SLNEW(DI),
     $         PopUp(DI), PopLow(DI),
     $         RadUpLow(DI), RadLowUp(DI),
     $         ColUpLow(DI), ColLowUp(DI),
     $         R_BRACKET(DI)

      ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_HYD_NLTE_TRA


      SUBROUTINE PRINT_NLTETRAPOP(PopLow, PopUp, DEPLOW, DEPUP)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      REAL*8, DIMENSION(DPN), INTENT(IN) :: PopLow, PopUp, DEPLOW, DEPUP

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      FILE_UNIT = 1653

      CALL OPEN_TO_APPEND(FILE_UNIT, NTP_FILE)

      DO DI = 1, DPN; WRITE(FILE_UNIT, '(E15.7,3(2x,E15.7))') PopLow(DI), PopUp(DI), DEPLOW(DI), DEPUP(DI); ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_NLTETRAPOP


      SUBROUTINE PRINT_LTE_LEV(LEV_NAM, LEV_NUM, LEV_POP)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           LEV_NAM

      INTEGER, INTENT(IN) ::                LEV_NUM

      REAL*8, DIMENSION(DPN), INTENT(IN) :: LEV_POP

      CHARACTER(:), ALLOCATABLE ::          LEV_NAME, SEP_FILE_NAME

      INTEGER ::                            SEP_FILE_UNIT

      INTEGER ::                            DI
 
      LEV_NAME = RM_CHAR(RM_CHAR(LEV_NAM, ' '), '.')

      IF (LEV_NAME .EQ. 'HMINUS1')   LEV_NAME = 'HMINUS'
      IF (LEV_NAME .EQ. 'ELECTRONS') LEV_NAME = 'ELECTR'

      SEP_FILE_NAME = TRIM(ADJUSTL(LTE_DIR_1//LEV_NAME))

      SEP_FILE_UNIT = LEV_NUM * 13

      CALL OPEN_TO_APPEND(SEP_FILE_UNIT, SEP_FILE_NAME)

      WRITE(SEP_FILE_UNIT, '(1x,A,5x,A,2x,A,10x,A,/)') 'LEV', 'NUM', 'DI', 'POP'

      DO DI = 1, DPN; WRITE(SEP_FILE_UNIT, '(A6,2(2x,I3),3x,E15.7)') LEV_NAME, LEV_NUM, DI, LEV_POP(DI); ENDDO

      CLOSE(SEP_FILE_UNIT)

      END SUBROUTINE PRINT_LTE_LEV


      SUBROUTINE PRINT_LTE_TRA(LowLevel, UpLevel, LineInd,
     $                         PopLow, PopUp, TranInt,
     $                         RadUpLow, RadLowUp)

      USE FILE_OPERATIONS
      USE STRING_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER*10, INTENT(IN) ::           UpLevel, LowLevel

      INTEGER, INTENT(IN) ::                LineInd

      REAL*8, DIMENSION(DPN), INTENT(IN) :: PopLow, PopUp

      REAL*8, DIMENSION(DPN), INTENT(IN) :: TranInt

      REAL*8, DIMENSION(DPN), INTENT(IN) :: RadUpLow, RadLowUp

      CHARACTER(:), ALLOCATABLE ::          FILE_NAME, TRAN_NAME, LowLev, UpLev

      INTEGER ::                            FILE_UNIT

      INTEGER ::                            DI

      LowLev = RM_CHAR(RM_CHAR(LowLevel, ' '), '.')

      IF (LowLev .EQ. 'HMINUS1') LowLev = 'HMINUS'

      UpLev = RM_CHAR(RM_CHAR(UpLevel, ' '), '.')

      TRAN_NAME = LowLev//'_'//UpLev

      FILE_NAME = TRIM(ADJUSTL(LTE_DIR_2//TRAN_NAME))

      FILE_UNIT = LineInd * 165

      CALL OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

      WRITE(FILE_UNIT, '(3x,A,6x,A,3x,A,2x,A,10x,A,2(14x,A),2(13x,A),/)') 'LL',
     $                                                                    'UL',
     $                                                                    'IND',
     $                                                                    'DI',
     $                                                                    'PU',
     $                                                                    'PL',
     $                                                                    'INT',
     $                                                                    'RUL',
     $                                                                    'RLU'

      DO DI = 1, DPN

         WRITE(FILE_UNIT, '(2(A6,2x),2(I3,1x),5(1x,E15.7))')
     $         LowLev, UpLev, LineInd, DI,
     $         PopUp(DI), PopLow(DI), TranInt(DI),
     $         RadUpLow(DI), RadLowUp(DI)

      ENDDO

      CLOSE(FILE_UNIT)

      END SUBROUTINE PRINT_LTE_TRA


      SUBROUTINE UPDATE_ATM_MOD(ELEC_CONC)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      REAL*8, DIMENSION(DPN), INTENT(IN) :: ELEC_CONC

      REAL*8 :: H, TEMP, ELEC_CONC_OLD, HEAVY_ELEM_CONC, V_TURB

      INTEGER :: I

      OPEN(UNIT = 1935, FILE = 'FAL_VD')

      CALL RM_FILE(UPD_ATM_MOD_FILE, '-f'); CALL OPEN_TO_APPEND(1936, UPD_ATM_MOD_FILE)

      DO I = 1, DPN

         READ(1935, *) H, TEMP, ELEC_CONC_OLD, HEAVY_ELEM_CONC, V_TURB

         WRITE(1936, '(F10.5,2x,F12.5,1x,ES15.7,1x,ES15.7,1x,F10.5)') H, TEMP, ELEC_CONC(I), HEAVY_ELEM_CONC, V_TURB

      ENDDO

      CLOSE(1935); CLOSE(1936)

      END SUBROUTINE UPDATE_ATM_MOD

      END MODULE LINPOP_AUX
