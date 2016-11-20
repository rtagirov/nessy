      module MOD_WRITMSC
      contains
      SUBROUTINE WRITMSC (ifl,string,n,char,iflag,idummy,ierr)

      implicit real*8(a-h,o-z)

      character string*N
	character*10 char

      ierr=1
      if (iflag.eq.-1) then
	   write (ifl,'(A10)') char
	   write (ifl,*) string
	   ierr=0
	else
	   write(6,*) ' MODFILE write mode unknown'
	   pause
	   stop
	endif

	return
	end subroutine
      end module