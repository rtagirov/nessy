C**********  MODULNAME: XRUDI     ******* 24/03/87  22.16.54.******    32 KARTEN
      MODULE MOD_XRUDI
      contains
      SUBROUTINE XRUDI (XJ,WAVENUM,XJC,XLAMBDA,ND,NF,L)
      use MOD_BNUE
      use MOD_TRADFUN
      use MOD_ERROR
C***  INTERPOLATION OF THE CONTINUUM RADIATION FIELD AT WAVENUM
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE

      implicit none
      integer,intent(in)  :: ND,NF,L
      real*8, intent(in)  :: WAVENUM, XLAMBDA(NF), XJC(ND,NF)
      real*8, intent(out) :: XJ
      real*8  ::             WLENG, A, B, H, P, TRAD, TRADA, TRADB
      integer ::             NA, NB, NH

      WLENG=1.d8/WAVENUM
      NA=1
      A=XLAMBDA(1)
      NB=NF
      B=XLAMBDA(NF)
      IF ((WLENG-A)*(WLENG-B) .GT. .0) THEN
         write (6,*) 'RUDIMENTAL TRANSITION OUTSIDE WAVELENGTH GRID'
         STOP 'ERROR'
         ENDIF
   10 IF ( NB-NA .EQ. 1) GOTO 12
      NH=(NA+NB)/2
      H=XLAMBDA(NH)
      IF ((WLENG-A)*(WLENG-H) .GT. .0) GOTO 13
      NB=NH
      B=H
      GOTO 10
   13 NA=NH
      A=H
      GOTO 10
   12 P=(WLENG-A)/(B-A)
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE
      TRADA=TRADFUN(XLAMBDA(NA),XJC(L,NA))
      TRADB=TRADFUN(XLAMBDA(NB),XJC(L,NB))
      TRAD=P*TRADB+(1.-P)*TRADA
!       if(TRAD==0) then
!         print '("xrudi: TRADAB = ",f0.4,X,f0.4)',TRADA,TRADB
!         print '("xrudi: NA,NB =  ",i0,x,i0)',NA,NB
!         print '("xrudi: XJC(L) = ",5(e10.3,X))',XJC(L,NA),XJC(L,NB)
!         print '("xrudi: TRAD =   ",f0.2)', TRAD
!       endif
      XJ=BNUE(WLENG,TRAD)

      RETURN
      END SUBROUTINE
      END MODULE
