      module MOD_DELPLA
      contains
C**********  MODULNAME: DELPLA    ******* 24/03/87  21.05.04.******     8 KARTEN
      FUNCTION DELPLA (T)
      use MOD_BNUE
C***  THIS FUNCTION IS USED BY SUBR. PRIFLUX FOR THE ESTIMATE OF TCOLOR
C***  AND CALCULATES THE RATIO B-NUE(XLAM1,T) /  B-NUE(XLAM2,T)

      implicit real*8(a-h,o-z)

      COMMON / COMDPL / XLAM1, XLAM2

      DELPLA=BNUE(XLAM1,T)/BNUE(XLAM2,T)

      RETURN
      END function
      end module
