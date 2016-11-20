      module MOD_ETLHIST
      contains
C**********  MODULNAME: ETLHIST   ******* 24/03/87  19.21.26.******   121 KARTEN
      SUBROUTINE ETLHIST (NLINE,LINE,LINEKEY,MODHIST,MAXHIST,JOBNUM,
     $      VDOP,VDOPOLD,LAST,LCARD)
C***  UPDATING THE MODEL HISTORY FOR PROGRAM ETL
      IMPLICIT NONE
      integer,intent(in   ) :: MODHIST(MAXHIST),MAXHIST,JOBNUM
      real*8, intent(in   ) :: VDOP,VDOPOLD
      integer,intent(  out) ::  NLINE,LAST
!      LOGICAL,intent(inout) :: LINEKEY(NLINE)
      LOGICAL,intent(inout) :: LINEKEY(:)
      character,intent(  out) :: LCARD*(*)
      character,intent(inout) :: line(:)*7
!      character :: STRING*(4*(nline+1)),WORD*8, CIND*3
      character :: STRING, WORD*8, CIND*3
   
      integer :: I,IA,IB,IPOS,IND,IFIRST,LASTIND,NL,NOLD
	lcard=' '
C***  CONVERSION OF LINE INDEX INTO I FORMAT, OMITTING ZERO'S (=RUDIMENTALS)
      NOLD=NLINE
      NLINE=0
      DO 30 NL=1,NOLD
c      DECODE (10,31,LINE(NL)) IND
	read (line(nl),31) cind
   31 FORMAT (4X,A3)
      IF (CIND .EQ. '000') GOTO 30
      NLINE=NLINE+1
      write(LINE(NLINE),'(A3)') CIND
      LINEKEY(NLINE)=LINEKEY(NL)
   30 CONTINUE
         
c      ENCODE (16,17,MODHIST(LAST+1)) JOBNUM
	write (lcard(1:16),17) jobnum
   17 FORMAT (1H/,I3,'. ETL',7X)
      LAST=16
     
      IF (VDOP.NE.VDOPOLD) THEN
c            ENCODE (16,3,MODHIST(LAST+1)) VDOP
	   write (lcard(17:32),3) vdop
    3    FORMAT (3X,5HVDOP= , F5.1,3X)
            LAST=LAST+16
      ENDIF
     
      IF (NLINE .EQ. 0) THEN
c         ENCODE (8,8,MODHIST(LAST+1))
         write (lcard(last+1:last+8),8)
    8    FORMAT (8HNO LINES)
         LAST=LAST+8
         GOTO 10
         ENDIF
     
C***  ENCODING THE TREATED LINE INDICES
      IA=1
      IB=1
   18 IF (IB.EQ.NLINE) GOTO 19
      IF (.NOT.(LINEKEY(IB).AND.LINEKEY(IB+1)) .AND. (LINEKEY(IB).OR.
     $     LINEKEY(IB+1))) GOTO 19
      IB=IB+1
      GOTO 18
     
   19 continue
      IF (LINEKEY(IA)) 
c         MODHIST(LAST)=8H   ETLA:
     $   write (lcard(last+1:last+8),'(A8)') '   ETLA:'
      IF (.NOT. LINEKEY(IA)) 
c      MODHIST(LAST)=8H FORMAL:
     $   write (lcard(last+1:last+8),'(A8)') ' FORMAL:'
	last=last+8
      IPOS=0
      STRING=' '
      DO I=IA,IB
        IF (I .GT. IA) THEN
          read (Line(i-1),4) LASTIND
          read (Line(i ),4) IND
    4     FORMAT (I3)
          IF (LASTIND .EQ. IND-1) THEN
            STRING(IPOS:IPOS)='-'
            IF (I .GT. IA+1) THEN
              IF (STRING(IPOS-4:IPOS-4) .EQ. '-') IPOS=IPOS-4
              ENDIF
            ENDIF
          ENDIF
        IPOS=IPOS+4
        write (word,'(A3)') line(i)
        STRING(IPOS-3:IPOS)=WORD
      ENDDO
      IFIRST=LAST+1
c      LAST=LAST+(IPOS-1)/8+1
      last=last+ipos
c      DO 2 I=IFIRST,LAST
c      K=8*(I-IFIRST)
c    2 ENCODE (8,9,MODHIST(I)) STRING(1+K:8+K)
      lcard(ifirst:last)=string(1:ipos)
c    9 FORMAT (A8)
      IB=IB+1
      IA=IB
      if (last.gt.100) last=100
      IF (IB.LE.NLINE) GOTO 18
     
   10 CONTINUE
c     IF (LAST.GT.MAXHIST) THEN
c           CALL REMARK ('NEW MODEL HISTORY TOO LONG')
c           STOP 'ERROR'
c           ENDIF
     
C***  SHORTENING OF THE PRESENT ETL HISTORY ENTRY, IF IT IS UNCHANGED
C***  SINCE THE LAST (FULL) ETL ENTRY
     
C***  FIRST: SEARCH FOR THE LAST FULL ETL ENTRY
C***  NOTE: IT IS RECOGNIZED ONLY FROM THE "L" AND SUBSEQUENT BLANKS IN THE
C***  SECOND WORD
c      DO 20 I=ISTART-2,1,-1
c      IF (MODHIST(I+1) .EQ. 8HL        ) THEN
cC***     I IS NOW THE START INDEX OF THE PRECEDING FULL ETL ENTRY
c         DO 23 J=ISTART+1,LAST
c         K=I+J-ISTART
c         IF (MODHIST(J) .NE. MODHIST(K)) GOTO 22
c   23    CONTINUE
c     
cC***     THE PRESENT ENTRY EQUALS THE OLD ONE
cC***     NOTE: THE BLANKS BEHIND "ETL" ARE NOW OVERWRITTEN
c         MODHIST(ISTART+1)=8HL    (UN
c         MODHIST(ISTART+2)=8HCHANGED)
c         LAST=ISTART+2
c         GOTO 22
c         ENDIF
c   20 CONTINUE
c   22 CONTINUE
c     
c      MODHIST(1)=LAST
c      CALL WRITMS (3,MODHIST,MAXHIST,7HMODHIST,-1,IDUMMY,IERR)

      RETURN
      END subroutine
      end module
