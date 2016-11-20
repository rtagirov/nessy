      MODULE MOD_HMINUS

      CONTAINS

      SUBROUTINE runSUB(CMD,JOB, timer, IDUMMY, reset)
      use MOD_WRSTART
      use MOD_STEAL
      use MOD_WRCONT
      use MOD_COMO
      use MOD_ETL
      character*(*) :: CMD
      character*7 :: JOB
      integer :: idummy, timer
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
      call finish(CMD,timer,IDUMMY,reset)
      END SUBROUTINE
      

      SUBROUTINE HMINUS_LTE

      use MOD_TICTOC
      use UTILS,only: cp
      use PARAMS_ARRAY

      USE FILE_OPERATIONS

      IMPLICIT NONE

      integer ipdim,    nbdim,  idummy,   itemp
      integer itsw
      integer lblaon,   inew,   ipmax,   nbmax
      integer timer,timer2
      integer nbinw, it
      real*8  almin, almax
      real*8  scaevt, absevt, scafac, absfac, scagri

      parameter (IPDIM=25,NBDIM=99)
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)
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
      IDUMMY=-1.  !***  TEST FOR TOO LARGE USE OF BLANK COMMON
      call TIC(timer)
      call TIC()
      rewind(99); write(99,'(A7)'),'wrstart'; JOB='WRSTART'
      PRINT *,' WIND STARTED: JOB=',JOB,'  tstart= ',0
      call tic(timer2)
      call runSUB('WRSTART',JOB, timer2, IDUMMY)    !***************** WRSTART ********************
      CALL runSUB('STEAL',JOB, timer2, IDUMMY)   !***************** STEAL **********************
      CALL runSUB('WRCONT',JOB, timer2, IDUMMY)  !********************* wrcont *****************
      itemp=itemp+1
      if (itemp.ge.2) itsw=1


   10 continue
      IF (JOB.EQ.'newline') INEW=1
      if (itsw.eq.1) itsw=0  !& CALL  TEMPCOR
      print *,'HMINUS: CYCLE STARTED'
      call tic(timer2)
      call runSUB('COMO', JOB, timer2, IDUMMY)
      call runSUB('ETL',  JOB, timer2, IDUMMY)
      print *,'HMINUS: TIME USED FOR CYCLE: '//writeTOC(timer2)
      IF (INEW.EQ.1) THEN
        print *,'HMINUS: TIME FOR CYCLE INCLUDING LINE BACKGROUND '//
     &              'RADIATION FIELD, JOB='//JOB 
        INEW=0
        ENDIF
      selectcase(JOB)
        case('wrcont','extrap') ; call runSUB('wrcont',JOB,timer2,IDUMMY)
        case('repeat')          ; goto 10;
      endselect
      PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
      REWIND 99; WRITE (99,'(A4)') 'exit'
      call BLANK_CHECK(IDUMMY)
      print *,'HMINUS: Elapsed time for total job: '//writeTOC(timer)

      CALL SYSTEM("echo '\n'END OF LTE RUN - $(date)'\n'")

      END SUBROUTINE


      SUBROUTINE FINISH(PROG,timer,IDUMMY,reset)
      USE MOD_TICTOC
      USE UTILS
      IMPLICIT NONE
      character :: PROG*(*)
      character :: P2*20
      integer   :: timer
      logical,optional :: reset
      integer,save :: counter = 1
      integer :: IDUMMY
      P2=PROG
      print *,'HMINUS-Finish: '//trim(P2)
     &  //'('//trim(adjustl(int2str(counter)))//')'
     &  //', time: '//trim(adjustl(writeTOC(timer)))
     &  //', total: '//trim(adjustl(writeTOC()))
      if(present(reset).and.reset) call tic(timer)
      counter=counter+1
      call BLANK_CHECK(IDUMMY)
      END SUBROUTINE
      !*** Check if the last common was changed.
      SUBROUTINE BLANK_CHECK(IDUMMY)
      IF(IDUMMY/=-1.)
     &   PRINT '("HMINUS: possible overflow of blank common ",i0)',IDUMMY
      END SUBROUTINE
      END MODULE
      

      PROGRAM HMINUS

      use MOD_WRSTART
      use MOD_STEAL
      use MOD_WRCONT
      use MOD_COMO
      use MOD_ETL
      use MOD_TICTOC
      use UTILS,only: cp
      use MOD_HMINUS
      use PARAMS_ARRAY

      USE COMMON_BLOCK
      USE FILE_OPERATIONS

      IMPLICIT NONE

      integer ipdim,    nbdim,  idummy,   itemp
      integer iadr9,    itsw
      integer lblaon,   inew,   modhist,  iadr,   ipmax,   nbmax
      integer timer,timer2
      integer nbinw
      real*8  almin, almax
      real*8  dumm1, dumm2, dumm3, dumm4, dumm5, dumm6, dumm7, dumm8
      real*8  dumm9, dumm10, dumm11, dumm12, dumm21
      real*8  scaevt, absevt, scafac, absfac, scagri

c                     1350           40500              30
      COMMON // DUMM1(NDIM,15),DUMM2(NDIM,NDIM,5),DUMM3(MAXATOM,6)
c                   644            25392                  16200
     $      ,DUMM4(NDIMP2,7),DUMM5(NDIMP2,NDIMP2,3),DUMM6(NDIM,NDIM,2)
c                  64000                60000
     $      ,DUMM7(NDDIM,NFDIM,4),DUMM8(NDDIM,MAXIND)
c                   22000           158              1200
     $      ,DUMM9(NFDIM,11),DUMM10(NFLDIM,2),DUMM11(NDDIM,15)
c                   36000
     $      ,DUMM12(NDDIM,NDIM,7)
c                    3000          3000          3000
     $      ,MODHIST(MAXHIST),IADR(MAXADR),IADR9(MAXADR)
c                   36000
     $      ,IDUMMY
c                              9750
      COMMON /COMIND/  DUMM21(MAXIND,13)
      parameter (IPDIM=25,NBDIM=99)
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW

      COMMON /LIBLFAC/ SCAFAC(NDDIM,NFDIM),ABSFAC(NDDIM,NFDIM)
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
      IDUMMY=-1.  !***  TEST FOR TOO LARGE USE OF BLANK COMMON

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
        CASE('lte_run')

             CALL HMINUS_LTE(); STOP

!             CALL SYSTEM("echo -e '\n'--------------------------------------------------------- LTE RUN
!     $ -----------------------------------------------------------------'\n'")

!             LTE_RUN = .TRUE.

!             REWIND(99); WRITE(99, '(A7)'), 'wrstart'

!             GOTO 100

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

      CALL WRSTART  !***************** WRSTART ********************
      call finish('WRSTART',timer2,IDUMMY,.true.)
!      call cp('RADIOC','RADIOC.wrstart')
!      call cp('RADIOL','RADIOL.wrstart')
!      call cp('RADIOCL','RADIOCL.wrstart')
!      call cp('POPNUM','POPNUM.wrstart')
      CALL STEAL (JOB)   !***************** STEAL **********************
      call finish('STEAL',timer2,IDUMMY)
!      call cp('RADIOC','RADIOC.steal')
!      call cp('RADIOL','RADIOL.steal')
!      call cp('RADIOCL','RADIOCL.steal')
!      call cp('POPNUM','POPNUM.steal')

      IF (JOB.NE.'wrcont') THEN
        PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
        REWIND 99
        WRITE (99,'(A4)') 'exit'
        call BLANK_CHECK(IDUMMY)
        print *,' Elapsed time for total wind-job:',toc(timer),' sec'
        STOP 'EXIT'
      ENDIF

    1 continue
      call tic(timer2)
      CALL WRCONT (JOB)  !********************* wrcont *****************
      itemp=itemp+1
      if (itemp.ge.2) itsw=1
      call finish('WRCONT',timer2,IDUMMY)

      IF (JOB.NE.'repeat'.AND.JOB.NE.'newline') THEN
         PRINT *,' NO NEW JOB TO BE ROUTED; JOB=',JOB
         REWIND 99
         WRITE (99,'(A4)') 'exit'
         IF (IDUMMY.ne.-1.) then
            PRINT *, ' possible overflow of blank common'
            print *,      IDUMMY
            ENDIF
            print *,'HMINUS: Elapsed time for total wind-job:',
     &              toc(timer),' sec'
         STOP 'EXIT'
      ENDIF

   10 continue
      IF (JOB.EQ.'newline') INEW=1
      if (itsw.eq.1) itsw=0  !& CALL  TEMPCOR
      print *,'HMINUS: CYCLE STARTED'
      call tic(timer2)

      CALL COMO        !**************COMO****************************
      call finish('COMO',timer2,IDUMMY)
      CALL ETL (JOB)   !**************ETL ****************************
!      if(first) call cp('RADIOC','RADIOC.etl')
!      if(first) call cp('RADIOL','RADIOL.etl')
!      if(first) call cp('RADIOCL','RADIOCL.etl')
!      if(first) call cp('POPNUM','POPNUM.etl')

      call finish('ETL',timer2,IDUMMY)
      CALL STEAL (JOB) !*************STEAL****************************
!      if(first) call cp('RADIOC','RADIOC.steal.2')
!      if(first) call cp('RADIOL','RADIOL.steal.2')
!      if(first) call cp('RADIOCL','RADIOCL.steal.2')
!      if(first) call cp('POPNUM','POPNUM.steal.2')
      first=.false.
      call finish('STEAL',timer2,IDUMMY)
      print *,'HMINUS: TIME USED FOR CYCLE: '//writeTOC(timer2)
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
      call finish('WRCONT',timer2,IDUMMY)

      PRINT *,'HMINUS: NO NEW JOB TO BE ROUTED; JOB=',JOB
      REWIND 99
      WRITE (99,'(A4)') 'exit'
      call BLANK_CHECK(IDUMMY)
      print *,'HMINUS: Elapsed time for total job: '//writeTOC(timer)

      CALL SYSTEM("echo '\n'END OF THE RUN - $(date)'\n'")

      END PROGRAM
