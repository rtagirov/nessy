      module string_operations

      implicit none

      contains

      function int2_to_char(intnum) result(charnum)

      integer, intent(in) :: intnum

      character(len = 2)  :: charnum

      write(charnum, '(i2)') intnum

      return

      end function int2_to_char


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
