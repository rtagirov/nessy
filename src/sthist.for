      MODULE MOD_STHIST

      CONTAINS

      SUBROUTINE STHIST(LINE,GAMMAL,GAMMAR,DELTAC,MODHEAD,JOBNUM,CORMAX,REDUCE,MODHOLD,tdiff)

!     UPDATING THE MODEL HISTORY FOR MAIN PROGRAM STEAL

      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)
     
      character LINE*120
      CHARACTER*104 MODHEAD, MODHOLD
      CHARACTER JOB*7, maxcor*8
      integer tdiff
     
      MAXCOR='UNDEF.  '
      IF (JOBNUM .GT. 1) THEN
         IF (CORMAX .GT. 1.d-100) THEN
            CORLOG=LOG10(CORMAX)
            write (MAXCOR,3) CORLOG
    3       FORMAT (F8.4)
         ENDIF
      ENDIF
      IF (GAMMAL .NE. .0) THEN
	   write (LINE,2) jobnum,
     $          DELTAC,GAMMAL,GAMMAR,MAXCOR
    2    FORMAT (1H/,I3,'. STEAL                DELTAC=',F5.1,
     $            ' GAMMAL=',F8.1,' GAMMAR=',F8.1,' COR.=',A8 )
c            LAST=LAST+11
      ELSE
         write (LINE,12) jobnum,maxcor
   12    FORMAT (1H/,I3,'. STEAL   COR.=',A8 )
c            LAST=LAST+4
      ENDIF
     
C***  SPECIAL ENTRY IF JOB IS STARTED FROM AN OLD MODEL
C      CALL JSYMGET (2LG7,JOB)
      rewind 99
      read (99,'(A7)') JOB
      IF ((JOB .EQ. 'wrstart') .AND. OLDSTART) THEN

	   write (LINE(31:100),10) MODHOLD(15:32)
   10    FORMAT ('   POPNUM FROM OLD MODEL (START: ',A18,')')
c         LAST=LAST+7
      ENDIF
     
      IF (REDUCE.NE.1.) THEN
         write (line(91:100),4) reduce
    4    FORMAT ('REDUCE=',F3.2)
      ENDIF

      write (line(101:120),6) tdiff
    6 FORMAT ('  RTime:',i8,' sec')

      RETURN
      END subroutine
      end module
