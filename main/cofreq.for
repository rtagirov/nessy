      module MOD_COFREQ
      contains
C**********  MODULNAME: COFREQ    ******* 24/03/87  20.39.20.******    32 KARTEN
      SUBROUTINE COFREQ (XR,XB,XMAX,ERXMIN,GDT,THIN)
      use MOD_FUNSCH
      use MOD_ERF_INF
      use MOD_REGULA
C***  CALCULATES THE CONFINING FREQUENCIES XR AND XB FOR ONE DIRECTION

      implicit real*8(a-h,o-z)

      COMMON / COMFUN / DELTAV,XMIN
      LOGICAL THIN
C***  EPS = ACCURACY DEMAND FOR REGULA
      DATA  EPS / 1.d-3 /
     
      IF (GDT .GE. 1.d0-2.d0*ERXMIN ) GOTO 3
      THIN=.FALSE.
C***  FOR INNER POINTS :
      ERFXR=GDT+ERXMIN
c      IF (ERFXR.le.0.) then
c         print *,' ERFXR,GDT,ERXMIN: ',ERFXR,GDT,ERXMIN
c         print *,' COFREQ: ERFXR, GDT set to EPS'
c         ERFXR=EPS
c         GDT=EPS
c      endif
      CALL REGULA (ERF_INF   ,XR,ERFXR,-XMAX,XMAX,EPS)
      XB=XMAX
     
      IF (XR+XMAX .LE. DELTAV) RETURN
C***  BOUNDARY ZONE
         DV2=DELTAV/2.d0
C***     OPTICAL THIN CASE
         IF (GDT .GE. 2.d0*(ERF_INF(DV2)-ERXMIN)-1.d0) GOTO 3
         XMIN=-XMAX
         CALL REGULA (FUNSCH,XR,GDT,XR-EPS,DV2,EPS)
         XB=DELTAV-XR
      RETURN
     
C***  OPTICAL THIN CASE
    3 THIN=.TRUE.
      XR=.0d0
      XB=.0d0
      RETURN
      end subroutine
      end module
