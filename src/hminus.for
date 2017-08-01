      MODULE MOD_HMINUS

      CONTAINS

      SUBROUTINE runSUB(CMD,JOB, timer, reset)
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
      

      SUBROUTINE HMINUS_LTE

      use MOD_TICTOC
      use UTILS,only: cp

      use file_operations

      IMPLICIT NONE

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

      CALL RM_FILE(NTP_FILE, '-vf')
      CALL RM_FILE(NTW_FILE, '-vf')

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

      END SUBROUTINE

      END MODULE
      

      PROGRAM HMINUS

      use MOD_WRSTART
      use MOD_STEAL
      use MOD_WRCONT
      use MOD_COMO
      use MOD_ETL
      use MOD_TICTOC
      use UTILS, only: cp
      use MOD_HMINUS

      use common_block
      use file_operations

      IMPLICIT NONE

      integer itemp
      integer itsw
      integer lblaon, inew, ipmax, nbmax
      integer timer,timer2
      integer nbinw
      real*8  almin, almax

      real*8 :: wrcont_start, wrcont_finish, wrcont_time
      real*8 :: como_start,   como_finish,   como_time
      real*8 :: etl_start,    etl_finish,    etl_time
      real*8 :: steal_start,  steal_finish,  steal_time
      real*8 :: cycle_start,  cycle_finish,  cycle_time

      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW

      character*7 :: JOB

      logical :: first=.true.

      LOGICAL :: VEL_FIELD_FILE_EXISTS

      LTE_RUN =             .FALSE.
      VEL_FIELD_FROM_FILE = .FALSE.

      INQUIRE(FILE = VEL_FIELD_FILE, EXIST = VEL_FIELD_FILE_EXISTS)

      IF (VEL_FIELD_FILE_EXISTS) THEN 

         VEL_FIELD_FROM_FILE = .TRUE.

         WRITE(*, '(A,1x,A,A,/)') 'Attention: the velocity field will be taken from the file', VEL_FIELD_FILE, '!'

      ENDIF

      itsw=0
      itemp=0
      ALMIN=0.    !***  default values for lineblanketing
      ALMAX=0.
      LBLAON=0

      call TIC(timer)
      call TIC()

  100 CONTINUE

      REWIND 99; READ (99,'(A7)') JOB

      WRITE(*, '(A,A)') 'JOB = ', JOB

      SELECTCASE(JOB)

        CASE('tmpcor');            GOTO 1
        CASE('wrcont');            GOTO 1
        CASE('repeat', 'newline'); GOTO 10
        CASE('wrstart');           CONTINUE
        CASE('lte')

             CALL HMINUS_LTE(); STOP

      CASE DEFAULT

          PRINT *,' OPTION NOT KNOWN; JOB=',JOB
          REWIND 99; WRITE (99,'(A4)') 'abort'
          PRINT *,' Elapsed time for total wind-job:',toc(timer),' sec'
          STOP 'ABORT'

      ENDSELECT

      CALL SYSTEM("echo '\n'START OF THE RUN - $(date)'\n'")

      CALL SYSTEM("echo DIRECTORY: $(pwd)'\n'")

      CALL tic(timer2)

      CALL RM_FILE('BROYDEN', '-vf')
      CALL RM_FILE('NEWBROYDEN', '-vf')

      CALL RM_FILE('etl_time.out', '-vf')

      CALL WRSTART; call finish('WRSTART', timer2, .true.)

      CALL STEAL(JOB); call finish('STEAL', timer2)

      IF (JOB.NE.'wrcont') THEN
        PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
        REWIND 99
        WRITE (99,'(A4)') 'exit'

        print *,' Elapsed time for total wind-job:',toc(timer),' sec'
        STOP 'EXIT'
      ENDIF

    1 continue

      call tic(timer2)

      call cpu_time(wrcont_start)

      CALL WRCONT(JOB)

      call cpu_time(wrcont_finish)

      wrcont_time = wrcont_finish - wrcont_start

      itemp = itemp + 1

      if (itemp .ge. 2) itsw = 1

      call finish('WRCONT',timer2)

      IF (JOB.NE.'repeat'.AND.JOB.NE.'newline') THEN
         PRINT *,' NO NEW JOB TO BE ROUTED; JOB=',JOB
         REWIND 99
         WRITE (99,'(A4)') 'exit'
         print *,'HMINUS: Elapsed time for total wind-job:', toc(timer), ' sec'
         STOP 'EXIT'
      ENDIF

   10 continue
      IF (JOB.EQ.'newline') INEW=1
      if (itsw.eq.1) itsw=0

      print *,'HMINUS: CYCLE STARTED'

      call cpu_time(cycle_start)

      call tic(timer2)

      call cpu_time(como_start)

      CALL COMO

      call cpu_time(como_finish)

      como_time = como_finish - como_start

      call finish('COMO', timer2)

      call cpu_time(etl_start)

      call system("echo -n $(date +%s) >> etl_time.out")

      CALL ETL(JOB)

      call system("echo ' '$(date +%s) >> etl_time.out")

      call cpu_time(etl_finish)

      etl_time = etl_finish - etl_start

      call finish('ETL', timer2)

      call cpu_time(steal_start)

      CALL STEAL (JOB)

      call cpu_time(steal_finish)

      steal_time = steal_finish - steal_start

      call cpu_time(cycle_finish)

      cycle_time = cycle_finish - cycle_start

      first = .false.

      call finish('STEAL', timer2)

      print*, 'HMINUS: TIME USED FOR CYCLE: '//writeTOC(timer2)

      if (lambda_iter == 0 .or. lambda_iter == 1) call system('rm -fv times.out')

      call open_to_append(231, 'times.out')

      write(231, '(I2,5(2x,F6.3))') lambda_iter, wrcont_time, como_time, etl_time, steal_time, cycle_time

      close(231)

      IF (INEW.EQ.1) THEN
        print *,'HMINUS: TIME FOR CYCLE INCLUDING LINE BACKGROUND '//
     &              'RADIATION FIELD, JOB='//JOB 
        INEW=0
        ENDIF
      selectcase(JOB)
        case('wrcont','extrap') ; goto 1;
        case('repeat')          ; goto 10;
      endselect
      
      call WRCONT(JOB)
      call finish('WRCONT',timer2)

      PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
      REWIND 99
      WRITE (99,'(A4)') 'exit'
      print *,'HMINUS: Elapsed time for total job: '//writeTOC(timer)

      CALL SYSTEM("echo '\n'END OF THE RUN - $(date)'\n'")

      END PROGRAM
