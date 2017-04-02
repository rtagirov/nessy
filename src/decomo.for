      module MOD_DECOMO

      contains

      SUBROUTINE DECOMO(LSOPA, LSINT, LBLANK)

!     DECODING INPUT OPTIONS FOR PROGRAM "COMO"

      implicit real*8(a-h,o-z)
     
      CHARACTER KARTE*80
      OPEN (1,FILE='CARDS')
      REWIND 1
      LSOPA=-1
      LSINT=-1
      LBLANK=0
     
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN

            PRINT *,' LINE BLANKETING OPTION DECODED COMO'
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:10) .EQ. 'PRINT INT ') THEN

            DECODE (80,7,KARTE) XL
    7       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN

            DECODE (80,7,KARTE) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            ENDIF
      GOTO 8
   66 CONTINUE
      CLOSE (1)
      RETURN

      END subroutine
      end module
