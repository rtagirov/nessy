      module MOD_ISRCHFGT

      contains

      function ISRCHFGT(NF, array, istep, val) result(idx)

      use MOD_ERROR

      real*8, dimension(NF), intent(in) :: array

      real*8 :: val, arr

      integer :: idx

      if (val .ge. array(NF)) then

         idx = NF + 1

         return

      endif

      i = 0

      do k = NF, 1, -istep

         arr = array(k)

         if (arr .gt. val) i = k

      enddo

      if (i .gt. 0) then

         idx = i

      else

         write (6,*) 'isrchfgt failed'

         pause

         call ERROR( 'isrchfgt failed' )

      endif

      return

      end function

      end module
