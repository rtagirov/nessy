      module MOD_PRIDAT

      contains

      SUBROUTINE PRIDAT(N,LEVEL,NCHARG, WEIGHT,ELEVEL,EION,EINST,
     $                  KODAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $                  NATOM,ELEMENT,NOM,ABXYZ,ATMASS)
C*******************************************************************************
C***  PRINTOUT OF THE ATOMIC DATA DECODED *************************************
C***  ADDITIONAL PRINTOUT OF RELATIVE ABUNDANCES (READ FROM INPUT CARDS)
C***  AND MASS FRACTIONS
C*******************************************************************************
      use MOD_OMEG
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in) :: NATOM
      DIMENSION NCHARG(N), WEIGHT(N),ELEVEL(N)
	real*4 wei
      DIMENSION EION(N)
      DIMENSION EINST(N,N)
      DIMENSION ALPHA(N),SEXPO(N)
      character*8 :: agaunt(N)
      DIMENSION NOM(N)
      DIMENSION ABXYZ(NATOM),ATMASS(NATOM)
c      DIMENSION ABXYZ(10),ATMASS(10)
      DIMENSION COCO(N,N,4)
      DIMENSION ALTESUM(4,N)
      DIMENSION KODAT(NATOM)
      DIMENSION UPSI(3)
      CHARACTER*3 KRUDI
      CHARACTER*4 KEYCOL(N,N)
      CHARACTER*10 LEVEL(N)
      CHARACTER*10 ELEMENT(NATOM)
     
C***  CALCULATION OF THE MEAN ATOMIC WEIGHT "ATMEAN"
      ATMEAN=0.

      DO 40 NA=1,NATOM
      ATMEAN=ATMEAN+ABXYZ(NA)*ATMASS(NA)
   40 CONTINUE
     


      IND=0
      PRINT 9
!    9 FORMAT(1H1,//,20X,'A T O M I C   D A T A   U S E D :',/,20X,33('='))
    9 FORMAT(/,20X,'A T O M I C   D A T A   U S E D :',/,20X,33('='))
     
C***  ELEMENTS ---------------------------------------------------------
      DO 29 NA=1,NATOM
      ABNA=ABXYZ(NA)
      ABLOG=LOG10(ABNA)
      FRACM=ABNA*ATMASS(NA)/ATMEAN
      PRINT 10, ABNA,ELEMENT(NA),ABLOG,FRACM
   10 FORMAT (//,20X,33('|'),10X,'RELATIVE ABUNDANCE (BY NUMBER):',2X,
     $           G10.3,/,
     $       20X,33('|'),41X,12('='),/,
     $       20X,'|',31X,'|',/,
     $       20X,'|',11X,A10,10X,'|',33X,'LOG(AB)=  ',F6.2,/,
     $       20X,'|',31X,'|',/,
     $       20X,33('|'),27X,'MASS FRACTION:',2X,G10.3,/,
     $       20X,33('|'),41X,12('='),/)
     
C***  LEVELS -----------------------------------------------------------
      PRINT 11
   11 FORMAT (//,10X,'1. LEVELS :',/,10X,11('-'),/,
     $ '  NR      NAME       WEIGHT     ENERGY    CHARGE   IONIZATION'
     $ ,' POT.',/,
     $ 30X,'(KAYSER)',15X,'(KAYSER)',/)
      DO 2 J=1,N
      IF (NOM(J) .NE. NA) GOTO 2
	wei=WEIGHT(J)
      NW=nint(wei)
      PRINT 3, J,LEVEL(J),NW,ELEVEL(J),NCHARG(J),EION(J)
    3 FORMAT (I3,3X,A10,I10,F12.2,I10,F12.2)
    2 CONTINUE
     
C***  LINE TRANSITIONS  ------------------------------------------------
      DO IT=1,3
         UPSI(IT)=10.D0**((IT-1)*0.61979D0+3.77815D0)
      ENDDO
      PRINT 4,(UPSI(IT)/1000.,IT=1,3)
      PRINT *
    4 FORMAT (//,10X,'2. LINE TRANSITIONS :',/,10X,20('-'),//,
     $ ' IND UP LOW         UP    LOW',5X,'VAC/AIR-WAVELENGTH    ',
     $ 'A(UP-LOW)',24X,'COLLISIONAL  TRANSITIONS ',/,
     $ '     (INDICES)   (CONFIGURATIONS)    (ANGSTROEM)   (PER SECOND)'
     $ ,11X,'(KEYWORD)',5X,'(UPSILON FOR ',3F6.0,' KK)',/)
C     $ ,11X,*(KEYWORD)*,10X,*(COEFFICIENTS)*,/)
      DO 1 I=2,N
      IF (NOM(I) .NE. NA) GOTO 1
      IM=I-1
      DO 17 J=1,IM
      IF ((NOM(J) .NE. NOM(I)) .OR. (NCHARG(J) .NE. NCHARG(I))) GOTO 17
      IND=IND+1
      KRUDI='   '
      IF (EINST(J,I) .EQ. -2.) KRUDI='RUD'
      WLENG=1.E8/(ELEVEL(I)-ELEVEL(J))
      IF (WLENG.GT.2000.) THEN
C***  LANG P.204
          ANU2=(ELEVEL(I)-ELEVEL(J))*(ELEVEL(I)-ELEVEL(J))
          ANAIR=6432.8D-8+2949810./(146.D+8-ANU2)+25540./(41.D+8-ANU2)
          WLENG=WLENG/(1.+ANAIR)
          ENDIF
      DO IT = 1,3
         TL = 10.**((IT-1)*0.61979+3.77815)
         ENE=1.
         CALL OMEG (N,TL,NCHARG,ELEVEL,EINST, OMEGA, I, J,
     $              EION,COCO,KEYCOL,WEIGHT,ALTESUM,NATOM,NOM,KODAT)

         UPSI(IT) = OMEGA/8.63D-6*WEIGHT(I)*SQRT(TL)
      ENDDO
c      IF (KEYCOL(I,J).EQ.'    ' .OR. KEYCOL(I,J).EQ.'NONE' .OR.
c     $    KEYCOL(I,J).EQ.'JEFF' .OR. KEYCOL(I,J)(:3).EQ.'UPS') THEN
c              PRINT 12,IND,I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),KRUDI,
c     $         KEYCOL(I,J)
c   12 FORMAT(3I4,3X,A10,3X,A10,F10.2,1P,E15.4,5X,A3,5X,A4)
c      ELSE
c              PRINT 5,IND,I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),KRUDI,
c     $         KEYCOL(I,J),(COCO(I,J,M),M=1,4)
c    5 FORMAT(3I4,3X,A10,3X,A10,F10.2,1P,E15.4,5X,A3,5X,A4,4(2X,E8.1))
c      ENDIF
      IF ((UPSI(1).GT.0.1).AND.(UPSI(2).GT.0.1).AND.(UPSI(3).GT.0.1))
     $   THEN
      PRINT 5,IND,I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),KRUDI,
     $         KEYCOL(I,J),(UPSI(IT),IT=1,3)
    5 FORMAT(3I4,3X,A10,3X,A10,F16.2,1P,E12.4,5X,A3,5X,A4,0P,3(X,F14.2))
      ELSE
      PRINT 12,IND,I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),KRUDI,
     $         KEYCOL(I,J),(UPSI(IT),IT=1,3)
   12 FORMAT(3I4,3X,A10,3X,A10,F16.2,1P,E12.4,5X,A3,5X,A4,3(X,E14.3))
      ENDIF

   17 CONTINUE
    1 CONTINUE
     
C***  CONTINUA  --------------------------------------------------------
      PRINT 7
    7 FORMAT(//,10X,'3. CONTINUA :',/,10X,13('-'),//,
     $ ' LOW UP         LOW     UP VAC/AIR-THRESHOLD     ',
     $ 'PHOTO CROSS-SECTION    SEATON COEFF.         GAUNT FACTOR'/
     $ ' (INDICES)   (CONFIGURATIONS)    (ANGSTROEM)     ',
     $ 17H  (10**-18 CM**2) ,'      ALPHA      S           FORMULA'/)
      DO 16 J=2,N
      IF (NOM(J) .NE. NA) GOTO 16
      JM=J-1
      DO 6 I=1,JM
C***  LEVELS MUST BELONG TO THE SAME ELEMENT
      IF (NOM(J) .NE. NOM(I)) GOTO 6
C***  CHARGES MUST DIFFER BY 1
      IF (NCHARG(J).NE.NCHARG(I)+1 ) GOTO 6
C***  UPPER LEVEL MUST BE GROUND STATE OF THAT ION
      IF (NCHARG(J) .EQ. NCHARG(J-1)) GOTO 6
      WLENG=1.E8/(EION(I)-ELEVEL(I))
      IF (WLENG.GT.2000.) THEN
C***  CONVERSION FROM VACUUM TO AIR WAVELENGTH
C***  LANG P.204
          ANU2=(ELEVEL(I)-ELEVEL(J))*(ELEVEL(I)-ELEVEL(J))
          ANAIR=6432.8D-8+2949810./(146.D+8-ANU2)+25540./(41.D+8-ANU2)
          WLENG=WLENG/(1.+ANAIR)
          ENDIF
      IF (ALPHA(I) .EQ. .0D0) THEN
            PRINT 8,
     $        I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),SEXPO(I),AGAUNT(I)
    8       FORMAT
     $      (2I4,3X,A10,3X,A10,F10.2,F15.3,11X,'HYDROGENIC',F5.1,10X,A8)
            ELSE
            PRINT 13,
     $         I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),ALPHA(I),SEXPO(I),
     $         AGAUNT(I)
   13       FORMAT
     $      (2I4,3X,A10,3X,A10,F10.2,F15.3,F18.3,F8.1,10X,A8)
      ENDIF
    6 CONTINUE
   16 CONTINUE
     
C***  SUM OF TRANSITIONS TO LTE LEVELS  --------------------------------
      DO 21 I=1,N
      IF ((NOM(I) .EQ. NA) .AND. (ALTESUM(1,I) .GT. .0)) GOTO 22
   21 CONTINUE
      GOTO 29
   22 PRINT 23
   23 FORMAT (//,10X,'4. SUM OF TRANSITIONS TO LTE LEVELS',/,
     $ 10X,35('-'),//,
     $ '     LOW         LOW            UP          A-SUM',
     $      '   TEMPERATURE FUNCTION'/
     $ ' (INDEX)   (CONFIG.)   (SUM RANGE)   (PER SECOND)',
     $      '     COEFF.1   COEFF.2'/)
      DO 25 LOW=1,N
      IF ((NOM(LOW) .EQ. NA) .AND. (ALTESUM(1,LOW) .GT. .0)) PRINT 24,
     $  LOW,LEVEL(LOW),ALTESUM(4,LOW),ALTESUM(1,LOW),ALTESUM(2,LOW),
     $   ALTESUM(3,LOW)
   24 FORMAT (I8,2X,A10,6X,A8,1P,E15.3,0P,F12.4,F10.4)
   25 CONTINUE
   29 CONTINUE
     
      RETURN
      END subroutine
      end module
