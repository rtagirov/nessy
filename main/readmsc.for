      MODULE MOD_READMSC
      contains
      SUBROUTINE READMSC (ifl,carray,n,char,ierr)
!	@see READMS but for characters
      !implicit real*8(a-h,o-z)
      IMPLICIT NONE
      integer,intent(in) 	  :: ifl,n
      character*10,intent(in)   :: char
      integer,intent(out) 	  :: ierr
      character*(n),intent(out) :: carray
	
	character*10 		  :: cread,help

      ierr=1
      READ (ifl,'(A10)') cread
	do while (cread(1:1).eq.' ')
	   help=cread(2:10)
	   cread=help
	enddo
      if (cread.ne.char) then
	   write (6,*) 'READMSC: KEYWORD MISMATCH'
	   write (6,*) cread,char
c	   pause
	   stop
	endif
	READ (ifl,*) carray
	ierr=0

	return
	end SUBROUTINE
        END MODULE
