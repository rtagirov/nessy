      module utils

      implicit none

      contains

      subroutine error(msg, p)
!     general error routine to terminate the program after last message has been printed

      use ifcore

      implicit none

      logical,       intent(in), optional :: P

      character*(*), intent(in) ::           msg

      logical ::                             P_

      P_ = .false.

      if (present(P)) P_ = P

      print '("error: ",$)'

      print*, msg

      CALL TraceBackQQ(MSG, USER_EXIT_CODE= -1_4)

      if (P_) pause 'PAUSED'

      stop 'error: abnormal program termination'

      end subroutine error

      !****************************************************************
      !*** asserts that TEST is true. If not print MSG and exit program
      !*** If WARN is present and /= 0 print the message and WARN, do not exit
      !*** If the opional argument P is set to true, pause after the 
      !*** printout and before the call to error if WARN==0
      subroutine assert(TEST,MSG,WARN,P)
        logical,      intent(in) :: TEST
        integer,      intent(in),optional :: WARN
        logical,      intent(in),optional :: P
        character*(*),intent(in),optional :: MSG
        integer :: WARN__
        if(TEST) return
        WARN__ = 0
        if(present(WARN)) WARN__= WARN
        if(WARN__ >= 0) then
          print '("assert: Warning ",i0," : Reason: ",A)',WARN__,MSG
          if(present(P).and.P) pause
        else
          if(WARN__/=0)print '("assert: Error Level set to ",i0)',WARN__
          call ERROR('assert: failed: Reason: '//MSG,P)
        endif
      end subroutine assert

      pure function str2int(X)
        integer :: str2int
        character*(*),intent(in) :: X
        read (X,*), str2int
      end function str2int

      pure function str2real(X)
        real*8 :: str2real
        character*(*),intent(in) :: X
        read (X,*), str2real
      end function str2real

      pure function int2str(X)
        integer,intent(in):: X
        character*10 int2str
        write(int2str,'(i10)') X
      end function int2str

      !*** returns a free fileunit
      !*** not threadsafe
      function getFileUnit(start)
      implicit none;
      integer,optional :: start
      integer :: getFileUnit
      integer :: i
      logical :: lOPENED
      integer :: start_
      start_=10
      if(present(start)) start_=start
      do i=start_,1000
        INQUIRE(unit=i,opened=lOPENED)
        if(.not.lOPENED) then
          getFileUnit=i
          return
        endif
      enddo
      call error('utils:getFileUnit: no open unit found')
      end function getFileUnit
      
      !*** Copy a file from FROMFILE to TOFILE using the opional OPTS
      !*** and set the returnstatus to the opional status
      subroutine cp(FROMFILE,TOFILE,OPTS,status)
      character*(*),intent(in) :: FROMFILE,TOFILE
      character*(*),intent(in),optional :: OPTS
      integer,optional,intent(inout) :: status
      integer :: istat
      integer,external :: system
      istat = system('cp '//OPTS//' '//FROMFILE//' '//TOFILE)
      if(present(status)) then
        status=istat
      else
        if(istat/=0)
     &     print '("UTILS:CP: Warning: system returned ",i0)',istat
      endif

      end subroutine

      end module
