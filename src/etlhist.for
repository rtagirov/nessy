      module MOD_ETLHIST

      contains

      SUBROUTINE ETLHIST (NLINE,LINE,LINEKEY,JOBNUM,VDOP,VDOPOLD,LAST,LCARD)

C***  UPDATING THE MODEL HISTORY FOR PROGRAM ETL
      IMPLICIT NONE
      integer,intent(in) :: JOBNUM
      real*8, intent(in) :: VDOP, VDOPOLD
      integer,intent(  out) ::  NLINE,LAST

      LOGICAL,intent(inout) :: LINEKEY(:)
      character,intent(  out) :: LCARD*(*)
      character,intent(inout) :: line(:)*7

      character :: STRING, WORD*8, CIND*3
   
      integer :: I,IA,IB,IPOS,IND,IFIRST,LASTIND,NL,NOLD
	lcard=' '
C***  CONVERSION OF LINE INDEX INTO I FORMAT, OMITTING ZERO'S (=RUDIMENTALS)
      NOLD=NLINE
      NLINE=0
      DO 30 NL=1,NOLD

	read (line(nl),31) cind
   31 FORMAT (4X,A3)
      IF (CIND .EQ. '000') GOTO 30
      NLINE=NLINE+1
      write(LINE(NLINE),'(A3)') CIND
      LINEKEY(NLINE)=LINEKEY(NL)
   30 CONTINUE
         

	write (lcard(1:16),17) jobnum
   17 FORMAT (1H/,I3,'. ETL',7X)
      LAST=16
     
      IF (VDOP.NE.VDOPOLD) THEN

	   write (lcard(17:32),3) vdop
    3    FORMAT (3X,5HVDOP= , F5.1,3X)
            LAST=LAST+16
      ENDIF
     
      IF (NLINE .EQ. 0) THEN

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

     $   write (lcard(last+1:last+8),'(A8)') '   ETLA:'
      IF (.NOT. LINEKEY(IA)) 

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

      last=last+ipos



      lcard(ifirst:last)=string(1:ipos)

      IB=IB+1
      IA=IB
      if (last.gt.100) last=100
      IF (IB.LE.NLINE) GOTO 18
     
   10 CONTINUE

      RETURN

      END subroutine

      end module
