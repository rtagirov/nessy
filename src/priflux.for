      module MOD_PRIFLUX

      contains

      SUBROUTINE PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,JOBNUM,
     $                    FWEIGHT,MODHEAD,AKEY )

!     PRINTOUT OF THE EMERGENT CONTINUUM FLUX

      use MOD_TRADFUN
      use MOD_TCOLOR
      use MOD_PRICOLR
      use CONSTANTS,only:CLIGHT_CGS,PI

      use phys

      implicit real*8(a-h,o-z)

      DIMENSION XLAMBDA(NF),EMFLUX(NF),FWEIGHT(NF),AKEY(NF)
      CHARACTER MODHEAD*104, IFSUM*8, JUMP*8, ITCOL*8
	character*5  CEDGE

C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI  (ERG/CM**2/S/STERAD/KELVIN**4)
      DATA STEBOL / 1.804696d-5 /
C***  HC = H * C  ( ERG * ANGSTROEM )
      DATA HC / 1.9864837d-8/

      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A104,  20X,'JOB NO.',I5,
     $//,30X,'E M E R G E N T  A S T R O P H Y S I C A L  F L U X',/,30X
     $,51('='),/,34X,43H( FLUX = PI * ASTROPHYS.FLUX = 4 * PI * H ),//,
     $ ' FREQUENCY      LAMBDA    LOG F-NUE   LOG F-LAMBDA  ',
     $ '   T-RAD     T-COLOR    F-NUE /     LOG INTEGRATED F-NUE ',
     $ '  CONT. JUMP',
     $/,'   INDEX     (ANGSTROEM)  (ERG/CM+2) (ERG/CM+2/S/A) ',
     $ ' (KELVIN)    (KELVIN)    B(TEFF)    (ERG/CM2/SEC)    ',
     $ '    (MAGNITUDES)',/)
     
C***  CALCULATION OF EFFECTIVE TEMPERATURE
      TEFF=(TOTOUT /STEBOL)**.25
C***  CALCULATION OF LUMINOSITY
C***  THE NUMERICAL CONSTANT DENOTES THE LOG OF 4*PI*STEBOL/LSOLAR
      ALUMI=2.d0*LOG10(RSTAR*TEFF*TEFF)-36.730d0
     
      FSUM=.0d0
C***  LOOP OVER ALL FREQUENCY POINTS  ----------------------------------
!	OPEN (UNIT=1, FILE="flux2.out",FORM="FORMATTED") 
      DO 2 K=1,NF
      XLAM=XLAMBDA(K)
      FK=EMFLUX(K)
      W=1.d8/XLAM
      IF (FK .LE. .0) THEN
            FNUE=.0
            FLAMBDA=.0
            ELSE
            FNUE=LOG10(FK)
            FLAMBDA=FNUE+LOG10(CLIGHT_CGS*W/XLAM)
            ENDIF
      TRAD=TRADFUN(XLAM,FK)
C***  CALCULATION OF TCOLOR FROM THE SLOPE OF EMFLUX BETWEEN K-1 AND K
      ITCOL='        '
      IF (K.GT.1) THEN
         CALL TCOLOR (XLAMBDA(K-1),XLAM,EMFLUX(K-1),FK,ITCOL)
         ENDIF
      FOVERB=FK/BNUE(XLAM,TEFF)
C***  INTEGRATION OF F-NUE FROM ZERO TO CURRENT WAVELENGTH
      FSUM=FSUM+FK*FWEIGHT(K)
      IF (FSUM .GT. 1.d-99 .AND. FSUM .LT. 1.d+99) THEN
         FSUMLOG=LOG10(FSUM)
         write (IFSUM,5) FSUMLOG
    5    FORMAT (F8.3)
         ELSE
         FSUMLOG=-99.
         FSUM=0.
         IFSUM='    -INF'
         ENDIF
C***  CONTINUUM JUMPS AT EDGE FREQUENCIES
      JUMP='        '
c      READ (AKEY(K),99) NEDGE
      decode (5,99,akey(k)) AEDGE
   99 FORMAT (A5)
c	print 99,aedge
c      IF (NEDGE .EQ. 'EDGE+') THEN
c*** by writing out AKEY in ASCII a bit can get lost
c      if (AEDGE .eq. 5HEDGE+.or.AEDGE.eq.5HFDGE+) then
	write (CEDGE,99) AEDGE
      if (CEDGE(2:5) .eq. 'DGE+') then
         IF (EMFLUX(K-1) .LE. .0 .OR. EMFLUX(K) .LE. .0) THEN
            JUMP='   UNDEF'
            ELSE
            write (JUMP,8) LOG10(EMFLUX(K)/EMFLUX(K-1))*2.5
    8       FORMAT (F8.3)
            ENDIF
         ENDIF
c	write (1,*) xlam,fnue,flambda,TRAD,ITCOL,FOVERB,IFSUM,JUMP
!	write (1,*) xlam,fnue,flambda
    2 PRINT 3, K,XLAM,FNUE,FLAMBDA,TRAD,ITCOL,FOVERB,IFSUM,JUMP
    3 FORMAT (I7,F15.2,2F13.3,F12.0,4X,A8,G14.3,6X,A8,10X,A8)	
!	close (unit=1)
C***  ENDLOOP  ---------------------------------------------------------
     
C***  CALCULATION OF H-LYMAN PHOTONS
      SNNUE=.0
      SINUE=.0
      DO 6 K=1,NF
      EMFLW=EMFLUX(K)*FWEIGHT(K)
      ANNUE=EMFLW*XLAMBDA(K)/HC
      SNNUE=SNNUE+ANNUE
      SINUE=SINUE+EMFLW
      IF (XLAMBDA(K) .LT. 911.55 ) THEN
         SNALF=SNNUE
         KMAX=K
      ENDIF
    6 CONTINUE
    7 CONTINUE
      IF (SNNUE .GT. 1.d-99 .AND. SNNUE .LT. 1.d+99) THEN
      XLAMMAX=(XLAMBDA(KMAX)+XLAMBDA(KMAX+1))*.5
C*** TEMPERATURE DEFINITION OF KRIS DAVIDSON APJ 317,760
      T0=SINUE/SNNUE/2.701178033d0/1.38062259d-16
      PRINT 44,SINUE,SNNUE,T0
   44 FORMAT(/,10X,'ASTROPHYSICAL-FLUX-INTEGRAL:',1PE10.3,10X,
     $             'PHOTONS:',E10.3,10X,'T0:',0PF8.0)
      ENDIF
      IF (SNALF .GT. 1.d-99 .AND. SNALF .LT. 1.d+99) THEN
      SNALF=LOG10(4.*PI*PI*RSTAR*RSTAR*SNALF)
      ELSE
      SNALF=-99.
      ENDIF
     
      PRINT 4,TEFF,ALUMI,SNALF,XLAMMAX,TOTIN/TOTOUT
    4 FORMAT (/,10X,'EFFECTIVE TEMPERATURE:',F8.0,' KELVIN ',
     $            20X,'NUMBER OF H-LYMAN PHOTONS PER SECOND:',/,
     $          10X,'LOG OF LUMINOSITY (SOLAR UNITS):',F7.3,
     $            19X,'LOG(N)=',F6.2,'  (INTEGRATED UP TO',F7.2,
     $            ' ANGSTROEM)',/,
     $          10X,'RADIATIVE ENERGY INPUT/OUTPUT:',F8.3)
     
      CALL PRICOLR (NF,XLAMBDA,EMFLUX,RSTAR)
     
      RETURN
      END subroutine
      end module
