      module MOD_TRADFUN
      contains
C**********  MODULNAME: TRADFUN   ******* 24/03/87  22.14.04.******    12 KARTEN
      FUNCTION TRADFUN (XLAMBDA,XJ)
C***  RADIATION TEMPERAURE IN KELVIN, FROM XJ = J-NUE (CGS UNITS)
C***  AND XLAMBDA IN ANGSTROEM
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C

      implicit real*8(a-h,o-z)

      DATA C1,C2 / 1.4388d8, 3.9724d+8 /

      TRADFUN=.0d0
      IF (XJ.LE..0) RETURN
      TRADFUN=C1/XLAMBDA/LOG(C2/XLAMBDA/XLAMBDA/XLAMBDA/XJ +1.d0)

      RETURN
      END function
      end module
