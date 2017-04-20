      MODULE MOD_DECSTE

      CONTAINS

      SUBROUTINE DECSTE(LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,IFRRA,ITORA,LSEXPO,
     $                  IFLUX,IDAT,LEVELPL,N,IPLOTF,NEWWRC,
     $                  NGAMR,NGAML,AGAMR,AGAML,LINE,LASTIND,
     $                  TPLOT,Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB,LBLANK)

      use common_block
      use file_operations

!     DECODING INPUT OPTIONS FOR STEAL.FOR

      IMPLICIT REAL*8(A - H, O - Z)

      CHARACTER KARTE*80
      DIMENSION LEVELPL(N)
      INTEGER   NGAMR(10),NGAML(10)
      DIMENSION AGAMR(10),AGAML(10)
      LOGICAL   LINE(LASTIND),PROLIB,NOTEMP,TPLOT
C***  COMMON / COMPLOT / TRANSFERS THE MAXIMUM OF THE Y-AXIS TO SUBR. PLOTT
      COMMON / COMPLOT / YMAX
     
C***  NO PRINT OPTIONS SET (DEFAULT)
      TEFFE=-1.
      LSRAT=-1
      LSPOP=-1
      LSEXPO=-1
      IHIST=-1
      IFLUX=-1
      IDAT=-1
      LBLANK=0
C***  OLDSTART=.F.:  POPNUMBERS ARE CALCULATED IN SUBR. POPZERO (DEFAULT)
C***  OLDSTART=.T.:  POPNUMBERS FROM OLD MODEL FILE
      OLDSTART=.FALSE.
C***  JOBMAX EXCEEDED (DEFAULT)
      JOBMAX=-1
C***  DEMANDED CONVERGENCE ACCURACY (DEFAULT)
      EPSILON=.005
C*** NO REDUCTION OF CORRECTIONS (DEFAULT)
      REDUCE = 1.
C***  NO POPNUMBER PLOT (DEFAULT)
      LEVELPL(1)=0
      NPLOT=0
C***  NO FLUX PLOT (DEFAULT)
      IPLOTF=-1
      PROLIB=.FALSE.
C***  NUMBER OF REPEAT ITERATIONS BETWEEN TWO WRCONT-JOBS (DEFAULT)
      NEWWRC=11
C***  NO SCHARMER AMPLIFICATION OF CORRECTIONS
      DO 21 I=1,10
      NGAMR(I)=999999999
      NGAML(I)=999999999
      AGAMR(I)=0.
      AGAML(I)=0.
   21 CONTINUE
      NGAMR(1)=0
      NGAML(1)=0
      MGC=1
      MGR=1
      MGL=1
      IPRICC=0
      IPRILC=0
C***  LINES THAT HAVE NOT BEEN TREATED IN THE RADIATION TRANSFER CALCULATIONS,
C***  I.E. WHICH ARE NOT QUOTED IN A 'LINE' OPTION CARD, ARE NOT SCHARMER-
C*** AMPLIFIED.
      DO 32 IN=1, LASTIND
      LINE(IN)=.FALSE.
   32 CONTINUE
C***  DISTANCE MODULUS (DEFAULT VALUE)
      Y0=.0
C***  NO SUPPRESSION OF TEMPERATURE CORRECTIONS (DEFAULT)
c*** when using the temperature correction procedurte:
c    this option should be switched off within steal
      NOTEMP=.TRUE.
C*** NO TEMPERATURE PLOT (DEFAILT)
      TPLOT=.FALSE.
     
      OPEN (1,FILE='CARDS')
      REWIND 1
    8 READ (1,4,END=66) KARTE
    4 FORMAT (A)

      IF ( KARTE(:10) .EQ. 'LINE BLANK') THEN
C                           ==========
            IF (KARTE(17:18).EQ.'OF') THEN
               LBLANK=0
               PRINT *,' STEAL: LINE BLANKETING OPTION OFF'
            ELSE IF (KARTE(17:18).EQ.'ON') THEN
               LBLANK=1
               PRINT *,' STEAL: LINE BLANKETING OPTION ON'
            ELSE IF (KARTE(17:18).EQ.'PR') THEN
               LBLANK=2
               PRINT *,' STEAL: LINE BLANKETING OPTION ON (print on)'
            ELSE IF (KARTE(17:18).EQ.'IN') THEN
               LBLANK=-1
               PRINT *,' STEAL: LB TABLE INPUT'
            ELSE
               PRINT *,' NO VALID OPTION'
            ENDIF
            GOTO 8
            ENDIF

      IF (KARTE(:6) .EQ. 'NEWWRC' ) THEN
C                         ======
            DECODE (80,9,KARTE) XL
    9       FORMAT (7X,F10.0)
            NEWWRC=NINT(XL)
            GOTO 8
            ENDIF
      IF (KARTE(:9) .EQ. 'PLOT FLUX') THEN
C                         =========
            IPLOTF=1
            GOTO 8
            ENDIF
      IF ( KARTE(:6) .EQ. 'PROLIB' ) THEN
C                          ======
         DECODE(80,13,KARTE) TEFFE,GRAD,ALDMDT,VINF,BET,Y0
   13    FORMAT(10X,6F10.3)
            IPLOTF=1
            PROLIB=.TRUE.
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT RATE' ) THEN
                          !===========
            DECODE (80,1,KARTE) XL,XFRRA,XTORA

!0        1        2         3         4         5
!23456789 12345678 123456789 123456789 123456789 1234567890
!RINT RATE IIIIIIIII         IIIIIIIIII        IIIIIIIIII
    1       FORMAT (11X,F10.0,10X,F10.0,8X,F10.0)
            LSRAT=NINT(XL)
            IFRRA=NINT(XFRRA)
            ITORA=NINT(XTORA)
            IF (LSRAT.EQ.0) LSRAT=1
            IF (IFRRA.GT.ITORA) ITORA=IFRRA
            PRINT '("DECSTE: RATE =",I5,":",I5,":",I5)', IFRRA,LSRAT,ITORA
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT HIST' ) THEN
C                          ==========
            IHIST=1
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT FLUX' ) THEN
C                          ==========
            IFLUX=1
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT DATOM' ) THEN
C                          ===========
            IDAT=1
            GOTO 8
            ENDIF
      IF (KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                         =========
            DECODE (80,3,KARTE) XL
    3       FORMAT (9X,F10.0)
            LSPOP=NINT(XL)
            IF (LSPOP.EQ.0) LSPOP=1
            PRINT *,' PRINT POP ',LSPOP,' DECODED'
            GOTO 8
            ENDIF
      IF (KARTE(:10) .EQ. 'PRINT EXPO' ) THEN
C                          ==========
            DECODE (80,1,KARTE) XL
            LSEXPO=NINT(XL)
            IF (LSEXPO .EQ. 0) LSEXPO=1
            PRINT *,' PRINT EXPO ',LSEXPO,' DECODED'
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'JOBMAX=' ) THEN
C                         =======
            DECODE (80,11,KARTE) XL
   11       FORMAT (7X,F10.9)
            JOBMAX=NINT(XL)
!      IF (JOBMAX .GE. 1000) THEN
!            PRINT *,'JOBMAX .GT. 1000: INVALID OPTION'
!            STOP 'ERROR'
!            ENDIF
            GOTO 8
            ENDIF
      IF (KARTE(:7) .EQ. 'EPSILON' ) THEN
C                         =======
            DECODE (80,14,KARTE) EPSILON
   14       FORMAT (8X,F10.0)
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'REDUCE' ) THEN
C                         ======
            DECODE (80,17,KARTE) REDUCE
   17       FORMAT (7X,F10.0)
            IF (REDUCE.EQ..0) REDUCE=.5
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'PLOT POP') THEN
C                         ========
            IF (NPLOT+1 .EQ. N) GOTO 8
            DECODE (30,5,KARTE) L1,L2,L3
    5 FORMAT (9X,3(I2,1X))
            IF (L1 .LE. 0) GOTO 8
            NPLOT=NPLOT+1
            ENCODE (8,2,LEVELPL(NPLOT)) L1,L2,L3
    2       FORMAT (3I2)
            LEVELPL(NPLOT+1)=0
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT CCORE') THEN
C                          ===========
            IPRICC=1
            GOTO 8
            ENDIF
      IF (KARTE(:11) .EQ. 'PRINT LCORE') THEN
C                          ===========
            DECODE (80,7,KARTE) XLC
    7       FORMAT (11X,F10.0)
            IPRILC=NINT(XLC)
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'GAMMA=') THEN
C                         ======
            MGC=MGC+1
            MGR=MGR+1
            MGL=MGL+1
            IF (MGC .GT. 10 .OR. MGR .GT. 10 .OR. MGL .GT. 10) THEN
               write (6,*) ('TOO MANY GAMMA OPTIONS')
               STOP 'ERROR'
               ENDIF
            DECODE (80,12,KARTE) GAMMA, ANG
   12       FORMAT (6X,F11.0,9X,F10.0)
            NGAMR(MGR)=ANG
            NGAML(MGL)=ANG
            AGAMR(MGR)=GAMMA
            AGAML(MGL)=GAMMA
            GOTO 8
            ENDIF







    6       FORMAT (7X,F10.0,9X,F10.0)
      IF (KARTE(:6) .EQ. 'GAMMAR' ) THEN
C                         ======
            MGR=MGR+1
            IF (MGR .GT. 10) STOP 'ERROR'
            DECODE (80,6,KARTE) GAMMAR,ANGR
            NGAMR(MGR)=ANGR
            AGAMR(MGR)=GAMMAR
            GOTO 8
            ENDIF
      IF (KARTE(:6) .EQ. 'GAMMAL' ) THEN
C                         ======
            MGL=MGL+1
            IF (MGL .GT. 10) STOP 'ERROR'
            DECODE (80,6,KARTE) GAMMAL,ANGL
            NGAML(MGL)=ANGL
            AGAML(MGL)=GAMMAL
            GOTO 8
            ENDIF
      IF ( KARTE(:4) .EQ. 'LINE' ) THEN
C                          ====
            DECODE (80,41,KARTE) IND1
   41       FORMAT(4X,I3)

            IF (KARTE(9:10) .EQ. 'TO') THEN

               CALL SYSTEM('grep LINE'//' '//datom_nlte//' '//'> temp.out')

               IND2 = NUM_OF_LINES('temp.out')

               CALL SYSTEM('rm temp.out')

!                  DECODE (80,43,KARTE) IND2
!   43             FORMAT(15X,I3)

            ELSE

                  IND2 = IND1

            ENDIF

            IF (IND2 .GE. IND1) THEN
                  INC=1
                  ELSE
                  INC=-1
                  ENDIF
            DO 45 IND=IND1,IND2,INC
            IF (IND .GT. LASTIND) THEN
               PRINT *,' ERROR STOP IN DECSTE:  IND  .GT. LASTIND'
               write (6,*) (' IND  .GT. LASTIND')
               STOP 'ERROR'
               ENDIF
            LINE(IND)=.TRUE.
   45       CONTINUE
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'OLDSTART') THEN
C                         ========
            OLDSTART=.TRUE.
            GOTO 8
            ENDIF
      IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                         ========
         NOTEMP=.TRUE.
         GOTO 8
         ENDIF
      IF (KARTE(:6).EQ.'PLOT T') THEN
C                       ======
            DECODE (80,15,KARTE) YMAX
   15       FORMAT (22X,F6.0)
            TPLOT=.TRUE.
         GOTO 8
            ENDIF
     
      GOTO 8

   66 CONTINUE
      CLOSE (1)
      RETURN

      END subroutine
      end module
