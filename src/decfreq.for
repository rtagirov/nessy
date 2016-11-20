      module MOD_DECFREQ
      contains
C**********  MODULNAME: DECFREQ   ******* 24/03/87  20.49.53.******    52 KARTEN
      SUBROUTINE DECFREQ (XLAMBDA,NF,NFDIM,TREF)
C*******************************************************************************
C***  DECODE THE FREQUENCY GRID (WAVELENGTHS IN A) FROM TAPE6 = FGRID
!***  Also returns the 
C*******************************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      real*8,intent(out) :: XLAMBDA(NFDIM),TREF
      integer,intent(out):: NF
      integer,intent(in) :: NFDIM
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
      CHARACTER KARTE*80

C***  IF NO REFERENCE TEMPERATURE IS SPECIFIED, TEFF IS DEFAULT
      TREF=TEFF
     
      NF=0
C***  DECODING INPUT CARDS FROM TAPE 7 = FGRID
      OPEN (71,FILE='FGRID')
    1 READ (71,6,END=66) KARTE
    6 FORMAT(A)
      IF (KARTE(:1) .EQ.'*' ) THEN
            PRINT 2,KARTE
    2       FORMAT (1X,A)
            GOTO 1
            ENDIF
      IF (KARTE(:5) .EQ. 'TREF=' ) THEN
            DECODE (80,16,KARTE) TREF
   16       FORMAT(5X,F20.0)
            GOTO 1
            ENDIF
CMH  	IF (KARTE(:5) .EQ. '   ' ) THEN            
	IF (KARTE(:5) .EQ. '   ' ) THEN            
            GOTO 66
      ENDIF
c***  end of changes by Margit Haberreiter


      NF=NF+1
      IF(NF.GT.NFDIM) THEN
            WRITE(6,*) 'FREQUENCY DIMENSION NFDIM INSUFFICIENT'
            STOP 'ERROR'
            ENDIF

      DECODE (80,5,KARTE) XLAMBDA(NF)

    5 FORMAT(F10.0)
      GOTO 1

   66 CONTINUE
      IF(NF.EQ.0) THEN
            WRITE(6,*) 'NO FREQUENCY SCALE ENCOUNTERED'
		  STOP 'ERROR'
			
            ENDIF
      DO 21 K=2,NF
		IF((XLAMBDA(K-1)-XLAMBDA(K))*(XLAMBDA(1)-XLAMBDA(NF)).LE..0) THEN
            WRITE(6,*) 'WAVELENGTH SCALE OUT OF SEQUENCE'
		 STOP 'ERROR'
	
		ENDIF

   21 CONTINUE
	   CLOSE(71)

      RETURN

      END subroutine

      end module
