      module MOD_PREF_SYN

      contains

      subroutine pref_syn(KARTE,N,ELEVEL,LINE,INDLOW,INDNUP,LASTIND,
     $                    VDOP,FMAX,FMIN,XMAX,VMAX,VSIDU,esca_wd,
     $                    DXOBS,NFOBS,XLAM,FREMAX,
     $                    NF,FNUEC)

!     CALLED FROM FORMAL FOR DECODING LINE OPTION CARDS
!     CALCULATES LINE QUANTITIES FOR DETECTED LINE
!     IN CASE OF LINE OVERLAP:  PREPARES ALSO ALL QUANTITIES FOR BLENDING LINES

      use mod_lipo

      implicit real*8(a - h, o - z)
     
      DIMENSION ELEVEL(N)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      CHARACTER*7 KARBL(5)
      CHARACTER KARTE*80
     
      
      if (karte(1:4).eq.'LINE') then
!     branch for line card driven wavelength selection

!        DECODE THE DESIRED LINE
         DECODE (7,1,KARTE) LCODE,LINE
 1       FORMAT (A4,I3)
         PRINT*, 'ACHTUNG: LINE = ', LINE
         IF ((LINE .GT. LASTIND) .OR. (LINE .LT. 1) .OR.
     $        (LCODE .NE. 4HLINE)) THEN
            PRINT 8, KARTE(1:10)
 8          FORMAT (5X,'UNRECOGNIZED LINE OPTION CARD:',A10)
            LINE=0
            RETURN
         ENDIF
     
C***     CALCULATE LINE QUANTITIES
         LOW=INDLOW(LINE)
         NUP=INDNUP(LINE)
         XLAM=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
         print *,' line vac-wavelength: ',xlam

      else if (karte(1:4).eq.'INTE') then


c***  branch for wavelength interval
         if (xlam.le.0. .or. FMAX.le.0.) then
            print *,' something wrong!'
            print *,' xlam=',xlam,' FMAX=',Fmax
            stop 'error PREF_SYN'
            endif
         LINE=1
      else
         print *,' this mode is not defined'
         stop ' error PREFO_SYN'
      endif

c*** test for number of frequency points
      if (NFOBS.le.0) then
            print *,' something wrong!'
            print *,' NFOBS=',nfobs
            stop 'error PREF_SYN'
            endif

C***     DEFINING THE FREQUENCY BAND
      IF (FMAX .LE. 0.0) FMAX=VMAX+XMAX+VSIDU+esca_wd

      IF (FMIN .LE. 0.0) FMIN=FMAX
    
C***     DEFINING INCREMENT OF THE OBSERVER'S FRAME FREQUENCY
      PRINT*, 'FMAX FMIN:', FMAX, FMIN, NFOBS
      DXOBS=-(FMAX+FMIN)/FLOAT(NFOBS-1)
      FREMAX=FMAX

C***  EMERGENT FLUX (FNUEC) IS OBTAINED BY INTERPOLATION AT THE LINE FREQUENCY
 
      return

      end subroutine

      end module
