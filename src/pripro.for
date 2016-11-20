      module MOD_PRIPRO
      contains
      SUBROUTINE PRIPRO (XLAM,VDOP,NFOBS,PROFILE,XOBS0,DXOBS,JOBNUM,
     $     VSINI,REFCON,MODHEAD,DLAM,LSPRO,IFIRST,NPHI,LPSTA,LPEND,
     $          XN,XMAX,JFIRST,JLAST,P,W,PHEAD,PROLIB,
     $       AKEY,WEIGHTI,WEIGHTJ,LEVELI,LEVELJ,EINST,FNUEC,RSTAR2,VMAX
     $     ,EQWI)
c
c          revised March 1992 in Goddard
c
C***********************************************************************
C***  PRINTOUT OF THE EMERGENT LINE PROFILES ON FILE OUTPUT
C***********************************************************************
      use CONSTANTS,only:CLIGHT_SI
      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION PROFILE(NFOBS),DLAM(NFOBS),P(JFIRST)
      LOGICAL PROLIB
      CHARACTER MODHEAD*104,LEVELI*10,LEVELJ*10,PHEAD*28,AKEY*80
     
C***  C = VELOCITY OF LIGHT IN KM/SEC
      real*8,parameter :: C = CLIGHT_SI/1d3
      !DATA C/2.9979E+5/
C***  CFACS: ALL CONSTANT FACTORS USED IN THE CALCULATION OF THE ABSOLUTE
C***  LINE EMISSION (ABLIEM) IN UNITS OF THE SOLAR LUMINOSITY LSUN=3.82E33 5I7/S
C***  CFACS = 1.E8*PI*4.*PI*C/LSUN
      DATA CFACS/3.098E-14/
     
      IF (IFIRST.EQ.1.OR.LSPRO.GT.0) PRINT 4, MODHEAD,JOBNUM
    4 FORMAT (10H1$$$$$$$$$,  A  ,4X,'AFTER JOB NO.',I5)
      IF (IFIRST.EQ.1.OR.LSPRO.GT.0.OR.
     $(LPSTA.NE.1.OR.JFIRST.NE.1.OR.JLAST.EQ.1) )
     $PRINT 1,INT(XN),XMAX ,NFOBS,NPHI,LPSTA,LPEND,JFIRST,JLAST
    1 FORMAT(/' NUMERICAL PARAMETERS: XN=',I3,
     $      ' , XMAX=',F4.1,' , NFOBS=',I4,' , NPHI=',I4,
     $      ' , PHI-NR.',I4,' TO',I4,' , P-NR.',I4,' TO',I4 )
      IFIRST=0
     
      F=1.499E-16*XLAM*XLAM*EINST*WEIGHTJ/WEIGHTI

      PRINT 6,AKEY,XLAM,LEVELI,LEVELJ,EINST,F,VSINI,VDOP
    6 FORMAT (//,10X,A10,5X,'LAMBDA =',F10.2,' (ANGSTROEM)',
     $ 5X,'FROM LEVEL ',A10,' (LOW)  TO LEVEL ',A10,' (UP)',
     $ /,10X,103('-'),/,
     $ 10X,'  A(UP-LOW) =',1P,E12.4,0P,' (1/SEC)',7X,' F=',F6.3,12X,
     $ 'VSINI/(KM/S)=',F5.0,7X,'VDOP/(KM/S)=',F5.0,/ )

      IF (PROLIB) PRINT 13, PHEAD
   13 FORMAT (10X,'PROLIB: ',A28)

      IF (LSPRO.GT.0) PRINT 8
    8 FORMAT (/,
     $ 20X,'DELTA-NUE           DELTA LAMBDA          LAMBDA        ',
     $ '      RELATIVE         EQUIV.WIDTH',/,
     $ 20X,'(DOPPLER UNITS)     (ANGSTROEM)         (ANGSTROEM)     ',
     $ '        FLUX           (ANGSTROEM)',/)
     
      KHALF=0
      KMAX=0
      PMAX=0.
      PHALF=0.
      EQWI=.0
      ASYM=.0
      ABEQWI=0.
      EMEQWI=0.
      DLAM2=(DLAM(2)-DLAM(1))/2.
     
      X=XOBS0+DXOBS
      IF (LSPRO.GT.0)
     $PRINT 2,    X,DLAM(1),DLAM(1)+XLAM,PROFILE(1),EQWI
C***  LOOP OVER ALL FREQUENCY POINTS
      DO 3 K=2,NFOBS
      X=XOBS0+K*DXOBS
         EQWI=EQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
         ASYM=ASYM+(PROFILE(K-1)*DLAM(K-1)+PROFILE(K)*DLAM(K))*DLAM2
      IF (PROFILE(K).LT.1.)
     $                  ABEQWI=ABEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      IF (PROFILE(K).GE.1.)
     $                  EMEQWI=EMEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      IF (PROFILE(K).GT.PMAX) THEN
            PMAX=PROFILE(K)
            PHALF=PMAX/2.+0.5
            XPEAK=DLAM(K)+XLAM
            KMAX=K
            ENDIF
      IF (PROFILE(K).GT.PHALF) KHALF=K
      IF (LSPRO.GT.0)
     $PRINT 2,    X,DLAM(K),DLAM(K)+XLAM,PROFILE(K),EQWI
    2 FORMAT ( 8X,F20.2,F20.2,F20.2,F20.3,F20.3)
    3 CONTINUE
      IF (LSPRO.LE.0)  PRINT 10, EQWI
   10 FORMAT (/,10X, 'EQUIVALENT WIDTH:',21X,F20.3,' ANG' )
     
C***  CALCULATION OF THE ABSOLUTE LINE EMISSION
            ABLINE=-CFACS*EQWI*FNUEC*RSTAR2/XLAM/XLAM
          ABLIEM=-CFACS*EMEQWI*FNUEC*RSTAR2/XLAM/XLAM
          ABLIAB=-CFACS*ABEQWI*FNUEC*RSTAR2/XLAM/XLAM
      PRINT 12, ABLINE,ABLIEM,ABLIAB
   12 FORMAT (/,10X,'ABSOLUTE LINE EMISSION (SOLAR LUM.):',11X,F11.3,
     $',  EMISSION:',F10.3,',  ABSORPTION:',F10.3 )
     
      IF (KHALF.EQ.0.OR.KHALF.EQ.NFOBS) THEN
      PRINT 9
    9 FORMAT(/,20X,'*** ATTENTION: ABSORPTION PROFILE! ***',/)
      GOTO 11
      ENDIF
      IF (KMAX.GT.0) THEN
C***  find blue half intensity point
         DO 14 K=2,KMAX
            IF (PROFILE(K).lt.PHALF) KHALFB=K
   14    enddo
         XHALF=DLAM(KHALF)-(DLAM(KHALF)-DLAM(KHALF+1)) /
     $      (PROFILE(KHALF)-PROFILE(KHALF+1))*(PROFILE(KHALF)-PHALF)
         XHALFB=-DLAM(KHALFB)+(DLAM(KHALFB)-DLAM(KHALFB+1)) /
     $      (PROFILE(KHALFB)-PROFILE(KHALFB+1))*(PROFILE(KHALFB)-PHALF)
         FWHM=XHALF+XHALFB
         XKM=XHALF/XLAM*C
         XKMB=XHALFB/XLAM*C
         XKMFWHM=FWHM/XLAM*C
C***  CALCULATION OF DIMENSIONLESS MOMENTS (LINE PROFILE)
C***  SOURCE: CASTOR ET AL. 1981, MNRAS 194, 547
         FMOM=C/XLAM/VMAX
         W0=-FMOM*EQWI
         W1=FMOM*FMOM*ASYM
         PRINT 7,        PMAX,XPEAK,XHALF,XKM,XHALFB,XKMB,FWHM,XKMFWHM,
     $                   W0,W1
    7    FORMAT(/,
     $          10X,'PROFILE PEAK:',38X,F7.3,' AT ',F8.2,' ANG'/
     $          10X,'RED HALF WIDTH AT HALF MAXIMUM:',18X,F8.2,
     $              ' ANG ',F8.2,' KM/SEC'/
     $          10X,'BLUE HALF WIDTH AT HALF MAXIMUM:',17X,F8.2,
     $              ' ANG ',F8.2,' KM/SEC'/
     $          10X,'FULL WIDTH AT HALF MAXIMUM - FWHM:',15X,F8.2,
     $              ' ANG ',F8.2,' KM/SEC'/
     $          10X,'DIMENSIONLESS MOMENTS (CASTOR ET',
     $          ' AL. 1981):   W0 = ',F9.5,/,56X,'W1 = ',F9.5)
      ENDIF
     
      REFCON=REFCON/FNUEC
      PRINT 5,REFCON
    5 FORMAT (/,10X,'REF. CONT. REL. TO CONT. FLUX OF MODEL:',F8.3)
     
   11 CONTINUE
      IF (JFIRST .EQ. JLAST)
     $PRINT*,'          INTENSITY PROFILE AT P(',JFIRST,') = ',P(JFIRST)
     $      ,'    FLUX WEIGHT=',W
     
      RETURN
      END subroutine
      end module
