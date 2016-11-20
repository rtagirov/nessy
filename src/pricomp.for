      module MOD_PRICOMP
      contains
C**********  MODULNAME: PRICOMP   ******* 06/08/87  20.09.12.******    75 KARTEN
      SUBROUTINE PRICOMP (NDIM,EINST,N,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     $                   STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     $                   INDLOW,INDNUP)
C***********************************************************************
C***  PRINTOUT OF THE CHEMICAL COMPOSITION OF THE WR MODEL ATMOSPHERE
C***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION NCHARG(N),NOM(N)
      DIMENSION ABXYZ(NATOM)
      DIMENSION ATMASS(NATOM),STAGE(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      CHARACTER*10 ELEMENT(NATOM)
      CHARACTER*2 SYMBOL(NATOM)
     
      PRINT 11, SYMBOL(1),ABXYZ(1)
!   11 FORMAT (1H1,
   11 FORMAT (/,10X,'C H E M I C A L  C O M P O S I T I O N',/,10X,
     $        38('='),3/,1X,' REL. ABUNDANCES (BY NUMBER):',
     $        5X,A2,5X,1PE11.4)
      DO 19 NA=2,NATOM
      PRINT 12, SYMBOL(NA),ABXYZ(NA)
   12 FORMAT (35X,A2,5X,1PE11.4)
   19 CONTINUE
     
      NION=1
      DO 20 NLEV=2,N
      IF ((NOM(NLEV) .NE. NOM(NLEV-1)) .OR.
     $    (NCHARG(NLEV) .NE. NCHARG(NLEV-1))) NION=NION+1
   20 CONTINUE
      PRINT 21, NATOM,NION,N,LASTIND
   21 FORMAT (///,1X,' STATISTICS:',5X,I3,'  ELEMENTS',/,18X,I3,'  IONS'
     $        ,/,18X,I3,'  LEVELS',/,17X,I4,'  LINE TRANSITIONS',///)
     
      PRINT 31
   31 FORMAT (1X,' INDEX',3X,'ELEMENT',6X,'ATOMIC MASS',3X,'IONS',3X,
     $        'MAIN ION',/,1X,52('-'))
      DO 39 NA=1,NATOM
      NAION=1
      DO 38 NLEV=NFIRST(NA)+1,NLAST(NA)
      IF (NCHARG(NLEV) .NE. NCHARG(NLEV-1)) NAION=NAION+1
   38 CONTINUE
      PRINT 32,NA,ELEMENT(NA),ATMASS(NA),NAION,SYMBOL(NA),INT(STAGE(NA))
   32 FORMAT (3X,I2,5X,A10,5X,F6.2,7X,I2,6X,A2,1X,I2)
   39 CONTINUE
     
      PRINT 41
   41 FORMAT (///,1X,' ELEMENT',3X,'ION',3X,'CHARGE',3X,'LEVELS',3X,
     $        'LINES',3X,'(RUD.)',/,1X,50('-'))
      DO 49 NA=1,NATOM
      N1=NFIRST(NA)
   44 NCH1=NCHARG(N1)
      IF (NCH1 .EQ. NCHARG(NLAST(NA))) THEN
          NION=NLAST(NA)-N1+1
      ELSE
C          NION= ISRCHNE(NLAST(NA)-N1+1,NCHARG(N1+1),1,NCH1)
         DO NS=N1,NLAST(NA)-1
	      IF (NCHARG(NS).EQ.NCH1) NION=NS-N1+1
	   ENDDO
      ENDIF
      NLINE=0
      NRUD=0
      DO 45 IND=1,LASTIND
      LOW=INDLOW(IND)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) THEN
          NLINE=NLINE+1
          IF (EINST(LOW,INDNUP(IND)) .EQ. -2.) NRUD=NRUD+1
      ENDIF
   45 CONTINUE
      PRINT 42, SYMBOL(NA),NCH1+1,NCH1,NION,NLINE,NRUD
   42 FORMAT (4X,A2,6X,I2,6X,I2,7X,I2,6X,I3,5X,I3)
      N1=N1+NION
      IF (N1 .LE. NLAST(NA)) GOTO 44
   49 CONTINUE
     
      RETURN
      END subroutine
      end module
