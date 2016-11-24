      MODULE STRING_OPERATIONS

      IMPLICIT NONE

      CONTAINS

!      FUNCTION INT_TO_CHAR(IntNum) RESULT(CharNum)

!      INTEGER, INTENT(IN) ::       IntNum

!      CHARACTER(:), ALLOCATABLE :: CharNum

!      IF (IntNum .LE. 9) THEN

!         WRITE(CharNum, '(I1)') IntNum

!         CharNum = '00'//CharNum

!      ELSEIF (IntNum .GE. 10 .AND. IntNum .LE. 99) THEN

!         WRITE(CharNum, '(I2)') IntNum

!         CharNum = '0'//CharNum

!      ELSEIF ()

!         WRITE(CharNum, '(I3)') IntNum

!      ENDIF

!      WRITE(CharNum, '(I6)') IntNum

!      CharNum = TRIM(ADJUSTL(CharNum))

!      RETURN

!      END FUNCTION INT_TO_CHAR


      FUNCTION RM_CHAR(in_str, target_char) RESULT(out_str)

      CHARACTER(*), INTENT(IN) :: in_str

      CHARACTER :: target_char

      CHARACTER(:), ALLOCATABLE :: out_str

      CHARACTER :: ch

      INTEGER :: j

      out_str = ' '

      DO j = 1, LEN_TRIM(in_str)

         ch = in_str(j : j)

         IF (ch .NE. target_char) out_str = TRIM(out_str) // ch

      ENDDO

      out_str = TRIM(ADJUSTL(out_str))

      RETURN

      END FUNCTION RM_CHAR

      END MODULE
