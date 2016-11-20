!****|*****************************************************************
!* A Timer. Set the timer with call TIC(), TOC() returns the #s since  
!* last call to TIC().
!*
!* secs=TICTOC() reads out the timer and then resets it.
!* writeTOC() writes the time in an human readable form into an char*9
!****|*****************************************************************
      module MOD_TICTOC
      implicit none
        integer,private,save :: TICTOC__ = 0
        integer,private,external :: time
      contains

      ! reset the timer
      subroutine TIC(TIMER)
        integer,optional,intent(out) :: TIMER
        if(present(TIMER)) then
          TIMER=time()
        else
          TICTOC__=time()
        endif
      end subroutine

      ! return the seconds since last call to TIC()
      function TOC(timer)
        integer,optional,intent(in) :: TIMER
        integer :: TOC
        if(present(timer)) then
          TOC=time()-timer
        else
          TOC=time()-TICTOC__
        endif
      end function

      ! return the seconds since last call to TIC() and reset the timer
      function TICTOC(timer)
        integer,optional,intent(inout) :: timer
        integer TICTOC
        TICTOC=TOC(timer)
        call TIC(timer)
      end function

      ! write the TOCs into a char*9 string. If TOC > 1000h, ignore seconds
      function writeTOC(timer)
      integer,optional,intent(in) :: timer
      character*9 :: writeTOC
      integer :: foo
      foo=TOC(timer)
      if(foo>=1000*3600) then    ! print hh"h"mm"m", ignore seconds
        write (writeTOC,'(i5,"h",i2.2,"m")')
     $    foo/3600, mod(foo/60,60)
      elseif (foo>=3600) then    ! print hh:mm:ss
        write (writeTOC,'(i3.2,":",i2.2,":",i2.2)')
     $    foo/3600, mod(foo/60,60), mod(foo,60)
      else        ! print mm:ss
        write (writeTOC,'(XXX , X ,i2.2,":",i2.2)')
     $    foo/60, mod(foo,60)
      endif
      end function
      function writeTICTOC(timer)
      integer,optional,intent(inout) :: timer
      character*9 :: writeTICTOC
      writeTICTOC=writeTOC(timer)
      call TIC(timer)
      end function
      end module
