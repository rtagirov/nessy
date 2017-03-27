      MODULE FILE_OPERATIONS

      IMPLICIT NONE

      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_1 =        'NLTE/LEV/'
      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_2 =        'NLTE/TRA/'
      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_3 =        'NLTE/RAT/'
      CHARACTER (LEN = 8),  PARAMETER :: LTE_DIR_1 =         'LTE/LEV/'
      CHARACTER (LEN = 8),  PARAMETER :: LTE_DIR_2 =         'LTE/TRA/'

      CHARACTER (LEN = 8),  PARAMETER :: LTE_LINE_INT_FILE = 'LINE_INT'
      CHARACTER (LEN = 8),  PARAMETER :: LTE_CONT_INT_FILE = 'CONT_INT'
      CHARACTER (LEN = 7),  PARAMETER :: LTE_LEV_POP_FILE =  'LEV_POP'

      CHARACTER (LEN = 3),  PARAMETER :: NRRM_FILE_NAME =    'RNR' ! NRRM = Net Radiative Rate Matrix
      CHARACTER (LEN = 3),  PARAMETER :: NCRM_FILE_NAME =    'CNR' ! NCRM = Net Collision Rate Matrix
      CHARACTER (LEN = 3),  PARAMETER :: NTRM_FILE_NAME =    'TNR' ! NTRM = Net Total Rate Matrix

      CHARACTER (LEN = 7),  PARAMETER :: fal_mod_file     =  'ATM_MOD'
      CHARACTER (LEN = 10), PARAMETER :: upd_fal_mod_file =  'ATM_MOD.UPD'

      character (len = 5),  parameter :: kur_mod_file =      'TABLE'

      CHARACTER (LEN = 12), PARAMETER :: VEL_FIELD_FILE =    'vel_field.in'

      CHARACTER (LEN = 5),  PARAMETER :: CONV_DIR =          'CONV/'

      CHARACTER (LEN = 6),  PARAMETER :: EDDI_FILE =         'EDDIES'

      CHARACTER (LEN = 10), PARAMETER :: NTP_FILE =          'NLTETRAPOP'

      CHARACTER (LEN = 7),  PARAMETER :: NTW_FILE =          'NLTEWAV'

      PUBLIC

      CONTAINS

      SUBROUTINE CLEAN_DIR(DIR)

      CHARACTER (LEN = *), INTENT(IN) :: DIR

!      CHARACTER(:), ALLOCATABLE :: DIRECTORY

!      DIRECTORY = TRIM(ADJUSTL(DIR))

!      CALL SYSTEM('rm -vrf'//' '//DIRECTORY//'*')
      CALL SYSTEM('rm -vrf'//' '//DIR//'*')

      END SUBROUTINE


      SUBROUTINE MKDIR(DIR)

      CHARACTER (LEN = *), INTENT(IN) :: DIR

!      CHARACTER(:), ALLOCATABLE :: DIRECTORY

!      DIRECTORY = TRIM(ADJUSTL(DIR))

!      CALL SYSTEM('mkdir -vp'//' '//DIRECTORY)
      CALL SYSTEM('mkdir -vp'//' '//DIR)

      END SUBROUTINE


      SUBROUTINE RM_FILE(FILE_NAME, OPTIONS)

      CHARACTER (LEN = *), INTENT(IN) ::           FILE_NAME

      CHARACTER (LEN = *), INTENT(IN), OPTIONAL :: OPTIONS

!      CHARACTER(:), ALLOCATABLE ::                 FILENAME

!      FILENAME = TRIM(ADJUSTL(FILE_NAME))

      CALL SYSTEM('rm '//OPTIONS//' '//FILE_NAME)

      END SUBROUTINE


      SUBROUTINE OPEN_TO_APPEND(FILE_UNIT, FILE_NAME)

      CHARACTER (LEN = *), INTENT(IN) :: FILE_NAME

      INTEGER, INTENT(IN) ::             FILE_UNIT

      LOGICAL ::                         FILE_EXISTS

      SELECTCASE(FILE_NAME)

            CASE('line_feautrier_matrix.out'); OPEN(UNIT = FILE_UNIT, FILE = FILE_NAME, FORM = 'FORMATTED'); RETURN

            CASE('full_feautrier_matrix.out'); OPEN(UNIT = FILE_UNIT, FILE = FILE_NAME, FORM = 'FORMATTED'); RETURN

            CASE DEFAULT
                
            INQUIRE(FILE = FILE_NAME, EXIST = FILE_EXISTS)

            IF (.NOT. FILE_EXISTS) OPEN(UNIT = FILE_UNIT, FILE = FILE_NAME, STATUS = 'NEW', ACTION = 'WRITE')

            IF (FILE_EXISTS)       OPEN(UNIT = FILE_UNIT, FILE = FILE_NAME, ACTION = 'WRITE', STATUS = 'OLD', POSITION = 'APPEND')

      ENDSELECT

      END SUBROUTINE OPEN_TO_APPEND


      SUBROUTINE REPLACE_PATTERN(TARGET_FILE, PATTERN_1, PATTERN_2, LEN_PAT_1, LEN_PAT_2)

      CHARACTER (LEN = *), INTENT(IN) :: TARGET_FILE

      INTEGER, INTENT(IN) :: LEN_PAT_1
      INTEGER, INTENT(IN) :: LEN_PAT_2

      CHARACTER (LEN = LEN_PAT_1), INTENT(IN) :: PATTERN_1
      CHARACTER (LEN = LEN_PAT_2), INTENT(IN) :: PATTERN_2

      CHARACTER(:), ALLOCATABLE :: PAT_1, PAT_2

      PAT_1 = TRIM(ADJUSTL(PATTERN_1))

      IF (PATTERN_2 .NE. '      -      ') THEN

         PAT_2 = TRIM(ADJUSTL(PATTERN_2))

      ELSE

         PAT_2 = PATTERN_2

      ENDIF

      CALL SYSTEM("IFS='%'"); CALL SYSTEM("sed -i 's/"//PAT_1//"/"//PAT_2//"/g' "//TARGET_FILE); CALL SYSTEM("unset IFS")

      END SUBROUTINE REPLACE_PATTERN


      FUNCTION NUM_OF_LINES(PATH_TO_FILE) RESULT(LINE_NUM)

      CHARACTER (LEN = *), INTENT(IN) :: PATH_TO_FILE

      INTEGER :: FILE_UNIT, LINE_NUM, IO

      FILE_UNIT = 13745
     
      OPEN(UNIT = FILE_UNIT, FILE = PATH_TO_FILE, ACTION = 'READ')

      LINE_NUM = 0; IO = 0

      DO WHILE (IO .EQ. 0)

         READ(FILE_UNIT, *, IOSTAT = IO)

         IF (IO .NE. 0) EXIT

         LINE_NUM = LINE_NUM + 1

      ENDDO

      CLOSE(FILE_UNIT)

      RETURN

      END FUNCTION NUM_OF_LINES


!      FUNCTION LAST_LINE_OF_FILE(PATH_TO_FILE) RESULT(LINE)
!
!      CHARACTER (LEN = *), INTENT(IN) :: PATH_TO_FILE
!
!      CHARACTER (LEN = *) :: LINE
!
!      INTEGER :: FILE_UNIT, LINE_NUM, I
!
!      FILE_UNIT = 12745
!
!      LINE_NUM = LINE_NUM(PATH_TO_FILE)
!
!      OPEN(UNIT = FILE_UNIT, FILE = PATH_TO_FILE, ACTION = 'READ')
!
!      DO I = 1, LINE_NUM; READ(FILE_UNIT, *) LINE; ENDDO
!
!      END FUNCTION LAST_LINE_OF_FILE


      FUNCTION READ_ATM_MOD(FILENAME, COL_NUM) RESULT(ARRAY)

      CHARACTER (LEN = *), INTENT(IN) ::   FILENAME

      CHARACTER (LEN = 1), INTENT(IN) ::   COL_NUM

      REAL*8, DIMENSION(:), ALLOCATABLE :: ARRAY

      INTEGER ::                           LINE_NUM

      INTEGER ::                           I, FILE_UNIT

      REAL*8, DIMENSION(:), ALLOCATABLE :: H, TEMP, ELEC_CONC
      REAL*8, DIMENSION(:), ALLOCATABLE :: HEAVY_ELEM_CONC, V_TURB

      REAL*8, DIMENSION(:), ALLOCATABLE :: rho, pressure, vturb

      character (len = 10000) ::           header

      logical ::                           file_exists

      real*8 :: w

      inquire(file = filename, exist = file_exists)

      if (.not. file_exists) stop 'Atmosphere model file has not been found. Abort.'

      FILE_UNIT = 1834

      if (filename .eq. fal_mod_file) then

          LINE_NUM = NUM_OF_LINES(FILENAME)

          IF (ALLOCATED(ARRAY)) DEALLOCATE(ARRAY)

          ALLOCATE(ARRAY(LINE_NUM))
          ALLOCATE(H(LINE_NUM))
          ALLOCATE(TEMP(LINE_NUM))
          ALLOCATE(ELEC_CONC(LINE_NUM))
          ALLOCATE(HEAVY_ELEM_CONC(LINE_NUM))
          ALLOCATE(V_TURB(LINE_NUM))

          OPEN(UNIT = FILE_UNIT, FILE = FILENAME, ACTION = 'READ')

          DO I = 1, LINE_NUM

             READ(FILE_UNIT, *) H(I), TEMP(I), ELEC_CONC(I), HEAVY_ELEM_CONC(I), V_TURB(I)

          ENDDO

          CLOSE(FILE_UNIT)

          SELECTCASE(COL_NUM)

              CASE('1'); ARRAY(1 : LINE_NUM) = H(1 : LINE_NUM)
              CASE('2'); ARRAY(1 : LINE_NUM) = TEMP(1 : LINE_NUM)
              CASE('3'); ARRAY(1 : LINE_NUM) = ELEC_CONC(1 : LINE_NUM)
              CASE('4'); ARRAY(1 : LINE_NUM) = HEAVY_ELEM_CONC(1 : LINE_NUM)
              CASE('5'); ARRAY(1 : LINE_NUM) = V_TURB(1 : LINE_NUM)

              CASE DEFAULT; STOP 'FUNCTION READ_ATM_MOD: COL_NUM ARGUMENT IS NOT RECOGNIZED. ABORT.'

          ENDSELECT

          DEALLOCATE(H)
          DEALLOCATE(TEMP)
          DEALLOCATE(ELEC_CONC)
          DEALLOCATE(HEAVY_ELEM_CONC)
          DEALLOCATE(V_TURB)

      endif

      if (filename .eq. kur_mod_file) then

          LINE_NUM = NUM_OF_LINES(filename) - 1 ! minus one because of the header

          IF (ALLOCATED(ARRAY)) DEALLOCATE(ARRAY)

          ALLOCATE(ARRAY(LINE_NUM))
          ALLOCATE(rho(LINE_NUM))
          ALLOCATE(temp(LINE_NUM))
          ALLOCATE(pressure(LINE_NUM))
          ALLOCATE(elec_conc(LINE_NUM))
          ALLOCATE(vturb(LINE_NUM))

          OPEN(UNIT = FILE_UNIT, FILE = FILENAME, ACTION = 'READ')

          READ(FILE_UNIT, *) header

          DO I = 1, LINE_NUM

             READ(FILE_UNIT, *) rho(i), temp(i), pressure(i), elec_conc(i), w, w, vturb(i)

          ENDDO

          CLOSE(FILE_UNIT)

          SELECTCASE(COL_NUM)

              CASE('1'); ARRAY(1 : LINE_NUM) = rho(1 : LINE_NUM)
              CASE('2'); ARRAY(1 : LINE_NUM) = temp(1 : LINE_NUM)
              CASE('3'); ARRAY(1 : LINE_NUM) = pressure(1 : LINE_NUM)
              CASE('4'); ARRAY(1 : LINE_NUM) = elec_conc(1 : LINE_NUM)
              CASE('7'); ARRAY(1 : LINE_NUM) = vturb(1 : LINE_NUM)

              CASE DEFAULT; STOP 'FUNCTION READ_ATM_MOD: COL_NUM ARGUMENT IS NOT RECOGNIZED. ABORT.'

          ENDSELECT

          DEALLOCATE(rho)
          DEALLOCATE(TEMP)
          DEALLOCATE(elec_conc)
          DEALLOCATE(pressure)
          DEALLOCATE(vturb)

      endif

      RETURN

      END FUNCTION READ_ATM_MOD

      END MODULE FILE_OPERATIONS
