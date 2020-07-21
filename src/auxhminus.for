      module auxhminus

      contains

      subroutine runSUB(CMD,JOB, timer, reset)

      use MOD_WRSTART
      use MOD_STEAL
      use MOD_WRCONT
      use MOD_COMO
      use MOD_ETL

      character*(*) :: CMD
      character*7 :: JOB
      integer :: timer
      logical,optional:: reset
      print *,'HMINUS: Call '//CMD
      selectcase(CMD)
        case('WRSTART'); call WRSTART
        case('WRCONT');  call WRCONT (JOB)
        case('STEAL');   call STEAL  (JOB)
        case('ETL');     call ETL    (JOB)
        case('COMO');    call COMO
        case default
        print *,'HMINUS: Unknwon CMD:"'//CMD//'"'
        STOP 'HMINUS: Unknwon CMD'
      endselect
      call finish(CMD,timer,reset)
      END SUBROUTINE
      

      subroutine hminus_lte

      use mod_tictoc
      use file_operations

      implicit none

      integer itemp
      integer itsw
      integer lblaon,   inew,   ipmax,   nbmax
      integer timer,timer2
      integer nbinw, it
      real*8  almin, almax

      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW

      character*7 :: JOB

      CALL SYSTEM("echo '\n'START OF LTE RUN - $(date)'\n'")

      CALL SYSTEM("echo DIRECTORY: $(pwd)'\n'")

!      CALL RM_FILE(NTP_FILE, '-vf')
!      CALL RM_FILE(NTW_FILE, '-vf')

      itsw=0
      itemp=0
      ALMIN=0.    !***  default values for lineblanketing
      ALMAX=0.
      LBLAON=0

      call TIC(timer)
      call TIC()
      rewind(99); write(99,'(A7)'),'wrstart'; JOB='WRSTART'
      PRINT *,' WIND STARTED: JOB=',JOB,'  tstart= ',0
      call tic(timer2)
      call runSUB('WRSTART',JOB, timer2)
      CALL runSUB('STEAL',JOB, timer2)
      CALL runSUB('WRCONT',JOB, timer2)
      itemp=itemp+1
      if (itemp.ge.2) itsw=1


   10 continue
      IF (JOB.EQ.'newline') INEW=1
      if (itsw.eq.1) itsw=0
      print *,'HMINUS: CYCLE STARTED'
      call tic(timer2)
      call runSUB('COMO', JOB, timer2)
      call runSUB('ETL',  JOB, timer2)
      print *,'HMINUS: TIME USED FOR CYCLE: '//writeTOC(timer2)
      IF (INEW.EQ.1) THEN
        print *,'HMINUS: TIME FOR CYCLE INCLUDING LINE BACKGROUND '//
     &              'RADIATION FIELD, JOB='//JOB 
        INEW=0
        ENDIF
      selectcase(JOB)
        case('wrcont','extrap') ; call runSUB('wrcont',JOB,timer2)
        case('repeat')          ; goto 10;
      endselect
      PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
      REWIND 99; WRITE (99,'(A4)') 'exit'

      print *,'HMINUS: Elapsed time for total job: '//writeTOC(timer)

      CALL SYSTEM("echo '\n'END OF LTE RUN - $(date)'\n'")

      END SUBROUTINE


      SUBROUTINE FINISH(PROG,timer,reset)
      USE MOD_TICTOC
      USE UTILS
      IMPLICIT NONE
      character :: PROG*(*)
      character :: P2*20
      integer   :: timer
      logical,optional :: reset
      integer,save :: counter = 1
      P2=PROG
      print *,'HMINUS-Finish: '//trim(P2)
     &  //'('//trim(adjustl(int2str(counter)))//')'
     &  //', time: '//trim(adjustl(writeTOC(timer)))
     &  //', total: '//trim(adjustl(writeTOC()))
      if(present(reset).and.reset) call tic(timer)
      counter=counter+1

      end subroutine

      end module
