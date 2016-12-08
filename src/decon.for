      module MOD_DECON
      contains
      SUBROUTINE DECON(LSOPA,LSINT,IFLUX,JOBMAX,
     $                 LPRIH,LPHNU,LPRIV,TEFF,LBLANK)
C***********************************************************************
C***  DECODING INPUT OPTIONS, CALLED FROM WRCONT *******************************
C***********************************************************************
C234567890 234567890 234567890 234567890 234567890 234567890 234567890 234567890 
      IMPLICIT REAL*8(A-H,O-Z)
      
      LOGICAL NOTEMP
      CHARACTER KARTE*80
     
C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      LPRIH=-1
      LPHNU=-1
      LPRIV=-1
      IFLUX=-1
      JOBMAX=-1
      LBLANK=0
      TEFF=-1.0d0
      NOTEMP=.FALSE.
     
      OPEN (1,FILE='CARDS')
      REWIND 1
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)
     
      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
               PRINT *,' WRCONT: LINE BLANKETING OPTION OFF'
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
               PRINT *,' WRCONT: LINE BLANKETING OPTION ON'
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
               PRINT *,' WRCONT: LINE BLANKETING OPTION ON (print on)'
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
               PRINT *,' WRCONT: LB TABLE INPUT'
            ELSE
               PRINT *,' NO VALID OPTIN DECODED'
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:10) .EQ. 'PRINT FLUX') THEN
C                          ==========
            IFLUX=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT INT ') THEN
C                          ==========
            read (karte,7) XL
    7       FORMAT (10X,F10.0)
            LSINT=NINT(XL)
            IF (LSINT.EQ.0) LSINT=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN
C                          ==========
            read (karte,7) XL
            LSOPA=NINT(XL)
            IF (LSOPA.EQ.0) LSOPA=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT F(R)') THEN
C                          ==========
            read (karte,7) XL
            LPRIH=NINT(XL)
            IF (LPRIH.EQ.0) LPRIH=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT HNU ') THEN
C                          ==========
            read (karte,7) XL
            LPHNU=NINT(XL)
            IF (LPHNU.EQ.0) LPHNU=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT VLJ ') THEN
C                          ==========
            read (karte,7) XL
            LPRIV=NINT(XL)
            IF (LPRIV.EQ.0) LPRIV=1
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'JOBMAX=') THEN
C                         =======
            read (karte,7) XL
    3       FORMAT (7X,F10.9)
            JOBMAX=NINT(XL)
!      IF (JOBMAX .GE. 1000) THEN
!            PRINT *,'JOBMAX .GT. 1000: INVALID OPTION'
!            STOP 'ERROR'
!            ENDIF
            GOTO 8
            ENDIF
      IF (KARTE(:5) .EQ. 'TEFF=') THEN
C                         =====
            read (karte,5) TEFF
    5       FORMAT (5X,F10.0)
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                         ========
            NOTEMP=.TRUE.
            GOTO 8
            ENDIF
     
      GOTO 8

   66 CONTINUE
      CLOSE (1)
      RETURN
     
      END subroutine
      end module
