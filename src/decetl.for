      MODULE MOD_DECETL

      CONTAINS

      SUBROUTINE DECETL(LSOPA, LSINT, VDOP, LINE, NLINE, LINEKEY, MAXIND, LBLANK)

      use file_operations

!*******************************************************************************
!     DECODES INPUT CARDS FOR MAIN PROGRAM "ETL"
!*******************************************************************************

      implicit none
      integer,intent(in   ) :: MAXIND
      real*8, intent(  out) :: VDOP
      integer,intent(  out) :: LSOPA,LSINT,NLINE,LBLANK
      logical,intent(inout) :: LINEKEY(MAXIND)
      character*7,intent(inout) :: LINE(MAXIND)
      CHARACTER :: KARTE*80
      LOGICAL :: ETLKEY
      real*8  :: XL
      integer :: IND,IND1,IND2,INC
      
      OPEN (1,FILE='CARDS')
      LSOPA=-1
      LSINT=-1
      NLINE=0
      ETLKEY=.FALSE.
      VDOP=.0
      LBLANK=0
     
    6 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            PRINT *,' ETL: LB OPTION DECODED  -  not active in this rout
     $ine'
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
            ENDIF
            GOTO 6
            ENDIF

      IF ( KARTE(:10) .EQ. 'PRINT INTL') THEN
C                           ==========
            DECODE (80,8,KARTE) XL
    8       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            GOTO 6
            ENDIF

      IF ( KARTE(:10) .EQ. 'PRINT OPAL') THEN
C                           ==========
            DECODE (80,8,KARTE) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            GOTO 6
      ELSEIF (KARTE(:5) .EQ. 'VDOP=' ) THEN
C                         =====
            DECODE (80,11,KARTE) VDOP
   11       FORMAT (5X,F10.0)
            GOTO 6
      ELSEIF ( KARTE(:4) .EQ. 'LINE' ) THEN
C                          ====
        read (KARTE,'(4X,I3)') IND1

        IF (KARTE(9:10) .EQ. 'TO') THEN

!          READ (KARTE,'(15X,I5)') IND2

!RT        This is a dirty fix for the hardcoding of the number of lines in the CARDS file.
!RT        Now the number of lines is calculated automatically using the DATOM file.

           CALL SYSTEM('grep LINE'//' '//datom_nlte//' '//'> temp.out')

           IND2 = NUM_OF_LINES('temp.out')

           CALL SYSTEM('rm temp.out')

!           PRINT*, 'ACHTUNG: IND2 = ', IND2

!           IF (IND2 .GT. 150) READ (KARTE,'(15X,I5)') IND2

        ELSE

           IND2 = IND1

        ENDIF

        IF (IND2 .GE. IND1) THEN

          INC=1

        ELSE

          INC=-1

        ENDIF

        PRINT '("decetl:LINE ",i0," TO ",i0, ", inc=",i0)', IND1, IND2, INC

        DO IND = IND1, IND2, INC

          NLINE = NLINE + 1

          IF (NLINE .GT. MAXIND) THEN
            PRINT *,' ERROR STOP IN DECETL: NLINE .GT. MAXIND'
            ! CALL REMARK ('NLINE .GT. MAXIND')
            STOP 'ERROR'
          ENDIF
          write (LINE(NLINE),'(4HLINE,I3)') ind
          LINEKEY(NLINE)=ETLKEY
        ENDDO
        GOTO 6
      ENDIF
      IF ( KARTE(:6) .EQ. 'FORMAL' ) THEN
C                          ======
            ETLKEY=.FALSE.
            GOTO 6
            ENDIF
      IF ( KARTE(:3) .EQ. 'ETL' ) THEN
C                          ===
            ETLKEY=.TRUE.
            GOTO 6
            ENDIF
      GOTO 6
66    CONTINUE
            CLOSE (1)

            RETURN

      END SUBROUTINE

      END MODULE
