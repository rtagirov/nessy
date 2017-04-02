      module MOD_TCOLOR

      contains

      SUBROUTINE TCOLOR (XL1,XL2,EM1,EM2,ITCOL)

!     CALCULATION OF TCOLOR FROM THE SLOPE OF EMFLUX BETWEEN XLAM1 AND XLAM2

      use MOD_REGULA

      use phys

      implicit real*8(a-h,o-z)
      character itcol*8

      COMMON / COMDPL / XLAM1, XLAM2

!     FUNCTION DELPLA IS FORMAL PARAMETER IN SUBR. REGULA
!     TCOL1,TCOL2 = MIN AND MAX COLOUR TEMPERATURE TO BE DETERMINED
!     EPS = ACCURACY OF TCOLOR DETERMINATION
!     DATA TCOL1,TCOL2,EPS/1000.d0,500000.d0,50.d0/
!  	  TCOL2 INCREASED TO 1000000

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
      return

      end subroutine

      FUNCTION DELPLA(T)

      use phys

!     THIS FUNCTION IS USED BY SUBR. PRIFLUX FOR THE ESTIMATE OF TCOLOR
!     AND CALCULATES THE RATIO B-NUE(XLAM1,T) /  B-NUE(XLAM2,T)

      implicit real*8(a - h, o - z)

      COMMON / COMDPL / XLAM1, XLAM2

      DELPLA = BNUE(XLAM1, T) / BNUE(XLAM2, T)

      return

      end function

      end module
