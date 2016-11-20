!***************************************************************
!*** 2006-01-17  Micha Schoell: a general error function to terminate the program after
!***							  a last message has been printed
!
      MODULE MOD_ERROR
      contains
      SUBROUTINE ERROR(MSG,P)
      USE IFCORE
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: MSG
      logical,      intent(in),optional :: P
      logical :: P_
      P_ = .false.
      if(present(P)) P_ = P
      PRINT '("error: ",$)'
      PRINT *,MSG
      CALL TraceBackQQ(MSG,USER_EXIT_CODE= -1_4)
      if(P_) pause 'PAUSED'
      !call print_my_stack()
      !PAUSE
      STOP 'ERROR: Abnormal Program termination'
      END SUBROUTINE
!        SUBROUTINE PRINT_MY_STACK
!        interface
!          function printstack(fd)
!            bind(c) printstack
!            integer, value :: fd
!            integer printstack
!          end function printstack
!        end interface
!        integer n
!        n = printstack(1)
!        print *, "printstack returned ", n
!        end subroutine print_my_stack
        END MODULE
