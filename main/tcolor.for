      module MOD_TCOLOR
      contains
C**********  MODULNAME: TCOLOR    ******* 24/03/87  22.14.03.******    28 KARTEN
      SUBROUTINE TCOLOR (XL1,XL2,EM1,EM2,ITCOL)
C***********************************************************************
C***  CALCULATION OF TCOLOR FROM THE SLOPE OF EMFLUX BETWEEN XLAM1 AND XLAM2
C***********************************************************************
      use MOD_BNUE
      use MOD_DELPLA
      use MOD_REGULA
      implicit real*8(a-h,o-z)
      character itcol*8
C***  COMMON /COMDPL / IS USED BY FUNCTION DELPLA (RATIO OF PLANCK FUNCTIONS)
      COMMON / COMDPL / XLAM1, XLAM2
C*** FUNCTION DELPLA IS FORMAL PARAMETER IN SUBR. REGULA
C***  TCOL1,TCOL2 = MIN AND MAX COLOUR TEMPERATURE TO BE DETERMINED
C***  EPS = ACCURACY OF TCOLOR DETERMINATION
CMH      DATA TCOL1,TCOL2,EPS/1000.d0,500000.d0,50.d0/
CMH	TCOL2 INCREASED TO 1000000
      DATA TCOL1,TCOL2,EPS/1000.d0,1000000.d0,50.d0/
      XLAM1=XL1
      XLAM2=XL2
         ITCOL='   UNDEF'
         BMIN=BNUE(XLAM2,TCOL1)
         IF (BMIN.GT.0.) THEN
         DFK=EM1/EM2
            DFKMIN=BNUE(XLAM1,TCOL1)/BNUE(XLAM2,TCOL1)
            DFKMAX=BNUE(XLAM1,TCOL2)/BNUE(XLAM2,TCOL2)
            IF (DFK.LT.DFKMAX.AND.DFK.GT.DFKMIN) THEN
               CALL REGULA (DELPLA,TCOL,DFK,TCOL1,TCOL2,EPS)
               write (ITCOL,10) TCOL
   10          FORMAT (F8.0)
             ENDIF
          ENDIF
      RETURN
      END subroutine
      end module