      MODULE MOD_READMS

      logical,private,save:: READMS_WCHA_FIRST=.true.

      contains

      SUBROUTINE READMS(ifl,array,n,keywrd_pass,ierr)
!     Read from the file #ifl into the array array(1:N).
!     Stops the program if the header
!     is not equal keywrd.
!     ierr is always 0 at return (used to be an error indicator?)
!     @see READMSC
!
!     micha - 2007-02-01: Added support for old-style radioc files with no WCHA* Data.
!               if no WCHA* data exist the array is filled with 0's

      use utils
      use MOD_FORMATS

      implicit none

      logical, save :: keywrd_read = .false.
      real*8 array(n)
      character*10, save:: cread
      character*20 keywrd
      character*(*) keywrd_pass
      integer :: i
      integer IFL,N,IERR,IOSTATUS

      ierr=1;
      keywrd=trim(adjustl(keywrd_pass))
      if(.not. keywrd_read) then
11      READ (ifl,'(A10)',iostat=IOSTATUS) cread
        if(IOSTATUS == -1) then  ! EOF Reached?
          write(6,*) 'EOF reached: ', keywrd
          if(keywrd(1:4) == 'WCHA') then
            array=0d0;
            return
          else
            call ERROR('READMS: EOF Reached unexpected!')
          endif
        endif
        keywrd_read=.true.
      endif
      if (.not.equalWithoutBlanks(cread,keywrd)) then
       if(cread(5:7)=='***'.and.cread(1:4)==keywrd(1:4)) then
        if(READMS_WCHA_FIRST) then
           write(6,'("readms: cread contains *s, ignore for now,'//
     $        ' cread=""",A," keywrd=""",A,"""")'), cread,keywrd
        endif 
       else
          if(keywrd(1:4).eq.'WCHA') then
            array(:)=0d0
            if(READMS_WCHA_FIRST) then
              write(6,*) 'READMS: WCHA not found - set to 0'  ! Warning: used to support old style
              READMS_WCHA_FIRST=.false.
            endif
            return
          else
            write (6,*) 'READMS: KEYWORD MISMATCH: "', keywrd,
     $                              '" != cread: "',cread,'"'

            call ERROR('READMS: KEYWORD MISMATCH')
            stop
          endif
        endif
	endif

	READ (ifl,*) (array(i),i=1,N)

      keywrd_read = .false.
	ierr=0


        return
	end SUBROUTINE

      SUBROUTINE READMS1 (ifl,scalar,keywrd,ierr)
!     Read from the file #ifl into the array array(1:N).
!     Stops the program if the header
!     is not equal keywrd.
!     ierr is always 0 at return (used to be an error indicator?)
!     @see READMSC
!
!     micha - 2007-02-01: Added support for old-style radioc files with no WCHA* Data.
!               if no WCHA* data exist the array is filled with 0's

      use utils

      use MOD_FORMATS

      implicit none
      real*8 scalar
      character*10, save:: cread
      character*(*) keywrd
      integer IFL,IERR,IOSTATUS

      ierr=1;
11    READ (ifl,'(A10)',iostat=IOSTATUS) cread
      if(IOSTATUS == -1) then  ! EOF Reached?
        write(6,*) 'EOF reached: ', keywrd
        call ERROR('READMS: EOF Reached unexpected!')
      endif

      if (.not.equalWithoutBlanks(cread,keywrd)) then
        write (6,*) 'READMS: KEYWORD MISMATCH: "', keywrd,
     $                            '" != cread: "',cread,'"'
        call ERROR('READMS: KEYWORD MISMATCH: "'
     $             //keywrd//'"!="'//cread//'"')
      endif
      READ (ifl,*) scalar
      ierr=0
      return
      end SUBROUTINE
      END MODULE
