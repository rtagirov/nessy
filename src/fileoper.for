      module file_operations

      implicit none

      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_1 =        'NLTE/LEV/'
      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_2 =        'NLTE/TRA/'
      CHARACTER (LEN = 9),  PARAMETER :: NLTE_DIR_3 =        'NLTE/RAT/'
      CHARACTER (LEN = 8),  PARAMETER :: LTE_DIR_1 =         'LTE/LEV/'
      CHARACTER (LEN = 8),  PARAMETER :: LTE_DIR_2 =         'LTE/TRA/'

      CHARACTER (LEN = 3),  PARAMETER :: NRRM_FILE_NAME =    'RNR' ! NRRM = Net Radiative Rate Matrix
      CHARACTER (LEN = 3),  PARAMETER :: NCRM_FILE_NAME =    'CNR' ! NCRM = Net Collision Rate Matrix
      CHARACTER (LEN = 3),  PARAMETER :: NTRM_FILE_NAME =    'TNR' ! NTRM = Net Total Rate Matrix

      CHARACTER (LEN = 7),  PARAMETER :: atm_mod_file     =  'atm.inp'

      CHARACTER (LEN = 12), PARAMETER :: VEL_FIELD_FILE =    'vel_field.in'

      CHARACTER (LEN = 5),  PARAMETER :: CONV_DIR =          'CONV/'

      CHARACTER (LEN = 6),  PARAMETER :: EDDI_FILE =         'EDDIES'

      CHARACTER (LEN = 10), PARAMETER :: NTP_FILE =          'NLTETRAPOP'

      CHARACTER (LEN = 7),  PARAMETER :: NTW_FILE =          'NLTEWAV'

      CHARACTER (LEN = 9),  PARAMETER :: atomic_data_file =  'datom.inp'

      public

      contains

      subroutine clean_dir(dir)

      CHARACTER (LEN = *), INTENT(IN) :: DIR

      CHARACTER(:), ALLOCATABLE :: DIRECTORY

      DIRECTORY = TRIM(ADJUSTL(DIR))

      CALL SYSTEM('rm -vrf'//' '//DIRECTORY//'*')

      end subroutine


      SUBROUTINE MKDIR(DIR)

      CHARACTER (LEN = *), INTENT(IN) :: DIR

      CHARACTER(:), ALLOCATABLE :: DIRECTORY

      DIRECTORY = TRIM(ADJUSTL(DIR))

      CALL SYSTEM('mkdir -vp'//' '//DIRECTORY)

      END SUBROUTINE


      SUBROUTINE RM_FILE(FILE_NAME, OPTIONS)

      CHARACTER (LEN = *), INTENT(IN) ::           FILE_NAME

      CHARACTER (LEN = *), INTENT(IN), OPTIONAL :: OPTIONS

      CHARACTER(:), ALLOCATABLE ::                 FILENAME

      FILENAME = TRIM(ADJUSTL(FILE_NAME))

      CALL SYSTEM('rm '//OPTIONS//' '//FILENAME)

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


      function num_of_lines(path_to_file) result(n)

      character (len = *), intent(in) :: path_to_file

      integer :: file_unit, n, io

      file_unit = 13745
     
      open(unit = file_unit, file = path_to_file, action = 'read')

      n = 0; io = 0

      do while (io == 0)

         read(file_unit, *, iostat = io)

         if (io .ne. 0) exit

         n = n + 1

      enddo

      close(file_unit)

      return

      end function

      function num_of_columns(path_to_file) result(n)

      character (len = *), intent(in) :: path_to_file

      integer :: fu, n

      call system('head -1'//' '//path_to_file//' '//'| wc -w > temp.out')

      fu = 13745

      open(unit = fu, file = 'temp.out', action = 'read')
      
      read(fu, *) n

      close(fu)

      call system('rm temp.out')
     
      return

      end function

      function read_atm_file_col(col) result(array)

      integer, intent(in) ::               col

      real*8, dimension(:), allocatable :: array

      integer ::                           nol

      integer ::                           i, file_unit

      real*8, dimension(:), allocatable :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10

      file_unit = 1834

      nol = num_of_lines(atm_mod_file)

      allocate(array(nol))

      allocate(c1(nol))
      allocate(c2(nol))
      allocate(c3(nol))
      allocate(c4(nol))
      allocate(c5(nol))
      allocate(c6(nol))
      allocate(c7(nol))
      allocate(c8(nol))
      allocate(c9(nol))
      allocate(c10(nol))

      open(unit = file_unit, file = atm_mod_file, action = 'read')

      selectcase(num_of_columns(atm_mod_file))

          case(4);  read(file_unit, *) (c1(i), c2(i), c3(i), c4(i),                      i = 1, nol) ! MURAM format

          case(5);  read(file_unit, *) (c1(i), c2(i), c3(i), c4(i), c5(i),               i = 1, nol) ! FAL format

          case(7);  read(file_unit, *) (c1(i), c2(i), c3(i), c4(i), c5(i), c6(i), c7(i), i = 1, nol) ! Kurucz 7  column format

          case(10); read(file_unit, *) (c1(i), c2(i), c3(i), c4(i), c5(i),
     $                                  c6(i), c7(i), c8(i), c9(i), c10(i), i = 1, nol)              ! Kurucz 10 column format

          case default; stop 'Function read_atm_file_col: Atmosphere model format is not recognized. Abort.'

      endselect

      close(file_unit)

      selectcase(col)

          case(1); array(1 : nol) = c1(1 : nol)
          case(2); array(1 : nol) = c2(1 : nol)
          case(3); array(1 : nol) = c3(1 : nol)
          case(4); array(1 : nol) = c4(1 : nol)
          case(5); array(1 : nol) = c5(1 : nol)
          case(6); array(1 : nol) = c6(1 : nol)
          case(7); array(1 : nol) = c7(1 : nol)

          case default; stop 'Function read_atm_file_col: col argument is not recognized. Abort.'

      endselect

      deallocate(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)

      return

      end function read_atm_file_col

      end module
