      MODULE MOD_READMSI
      contains
      SUBROUTINE READMSI (ifl,iarray,n,char,ierr)
      use MOD_ERROR
      use MOD_FORMATS
!      implicit real*8(a-h,o-z)
      implicit none
      integer :: IFL,N,IERR,I
      integer :: iarray(N)
	character*(*) char
      character*10 cread
      ierr=1
      READ (ifl,'(A10)') cread
      if (.not.equalWithoutBlanks(cread,char)) then
         write (6,*) 'READMSI: KEYWORD MISMATCH'
         call ERROR('READMSI: KEYWORD MISMATCH: "'//cread//
     $    '" != '//char//'"')
	endif
	READ (ifl,*) (iarray(i),i=1,N)
! D      if(cread == 'JOBNUM') print '("readmsi : JOBNUM=",i10)',iarray
	
      ierr=0
	return
	end SUBROUTINE
      SUBROUTINE READMSI1(ifl,iscalar,char,ierr)
      use MOD_ERROR
      use MOD_FORMATS
!      use IFCORE
!      implicit real*8(a-h,o-z)
      implicit none
      integer :: IFL,IERR
      integer :: iscalar
      character*10 cread
      character*(*) char

      ierr=1
      READ (ifl,'(A10)') cread
      if (.not.equalWithoutBlanks(cread,char)) then
         write (6,*) 'READMSI1: KEYWORD MISMATCH'
         call ERROR('READMSI1: KEYWORD MISMATCH: "'//cread//
     $    '" != '//char//'"')
      endif
      READ (ifl,*) iscalar
! D      if(cread == 'JOBNUM') then
! D        print '("readmsi1: JOBNUM=",i10)',iscalar
! D        if(iscalar==0)  CALL TraceBackQQ( USER_EXIT_CODE = -1)
! D     endif
      ierr=0
      return
      end SUBROUTINE

        end MODULE
