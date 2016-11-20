      module MOD_LINSCA
      contains
      SUBROUTINE LINSCA (XLAM,ELDEN,SCAFAC,ABSFAC)
ctest..
c         enhanced factor for HeII Ly-alpha (not active)
c
C  THIS SUBROUTINE INTERPOLATES THE RATIO OF LINE-SCATTERING TO
C  ELECTRON SCATTERING LINEARLY IN A TABLE
C
C  THE TABLE IS OUTPUT OF THE MONTE-CARLO PROGRAMM 
C     ==> IP AND DIMENSIONS HAVE TO BE UPDATED, IF THE TABLE IS CHANGED
      implicit real*8(a-h,o-z)
      real*8,intent(in):: xlam, elden
      real*8,intent(out):: scafac,absfac
      PARAMETER ( ONE = 1.D+0 )
      parameter (IPDIM=25,NBDIM=99)
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW

      IF (XLAM.LT.ALMIN .OR. LBLAON.eq.0) THEN
         SCAFAC=one
         ABSFAC=one
         RETURN
      ELSE IF (XLAM.LT.ALMAX) THEN
         Xbin=XLAM
         if (xbin.lt.911.7 .and. xbin.gt.900.) xbin=899.
         if (xbin.lt.227.9 .and. xbin.gt.220.) xbin=219.
         NLAM=Xbin/NBINW
         NLAM=NLAM+1
         IF (NLAM.GT.NBMAX) NLAM=NBMAX
      ELSE 
         SCAFAC=one
         ABSFAC=one
         RETURN
      ENDIF

      IF (NLAM.GT.NBDIM) STOP 'LINSCA'
      IF (NLAM.LT.1 .OR. NLAM.GT.NBMAX) STOP 'LINSCA'
      IF (ELDEN.GT.SCAGRI(1)) THEN
         SCAFAC=SCAEVT(1,NLAM)
         ABSFAC=ABSEVT(1,NLAM)
      ELSE IF (ELDEN.LT.SCAGRI(IPMAX)) THEN
         SCAFAC=SCAEVT(IPMAX,NLAM)
         ABSFAC=ABSEVT(IPMAX,NLAM)
      ELSE
         DO 1 I=1,IPMAX-1
            IF (ELDEN.LE.SCAGRI(I)) K1=I
 1       CONTINUE
         K2=K1+1
         DSCAGRI=SCAGRI(K1)-SCAGRI(K2)
         DSCAEVT=SCAEVT(K1,NLAM)-SCAEVT(K2,NLAM)
         SCAFAC=(ELDEN-SCAGRI(K1))/DSCAGRI*DSCAEVT+SCAEVT(K1,NLAM)
         DABSEVT=ABSEVT(K1,NLAM)-ABSEVT(K2,NLAM)
         ABSFAC=(ELDEN-SCAGRI(K1))/DSCAGRI*DABSEVT+ABSEVT(K1,NLAM)
ctest         if (xbin.lt.319. .and. xbin.gt.300.) scafac=scafac*10.
      ENDIF
      if (scafac.lt.one) stop ' scafac.lt.1 '
      if (absfac.lt.one) stop ' absfac.lt.1 '
      RETURN
      END subroutine
      end module
