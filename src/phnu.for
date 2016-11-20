      module MOD_PHNU
      contains
      SUBROUTINE PHNU (K,XL,LPHNU,ND,RADIUS,HNU,RSTAR)
      !*** Print the Eddington Flux as a function of Wavelength
      !** Does not change any variables.
      implicit none
      !IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in) :: K,LPHNU,ND
      real*8, intent(in) :: XL,RADIUS,HNU,RSTAR
      INTEGER :: I,NPRPT,L
      DIMENSION RADIUS(ND),HNU(ND)
      INTEGER LPT(100)
      IF (ND.LE.2) STOP 'ND.LE.2'
      IF (ND.GT.100) STOP 'ND.G.100'
      I=0
      IF (K.EQ.1) PRINT 10
   10 FORMAT (1H1/,
     $10X,'EDDINGTON FLUX AS A FUNCTION OF WAVELENGTH AND DEPTH',/,
     $10X,'===================================================='/)
CC      IF (K.EQ.1) PRINT *,RSTAR,' RC (CORE RADIUS IN CM)'
      NPRPT=(ND-2)/LPHNU
CC      IF (K.EQ.1) PRINT *,NPRPT,' PRINTED DEPTH POINTS'
      DO 111 L=2,ND-1
      IF(((L-1)/LPHNU)*LPHNU.NE.(L-1) .AND. L.NE.ND) GOTO 111
CC      IF (K.EQ.1) PRINT *,L,RADIUS(L)*RSTAR,RADIUS(L),' L, R, R/RC'
      I=I+1
      LPT(I)=L
  111 CONTINUE

      PRINT *,XL,' WAVELENGTH IN A'

C     DO 1 L=2,ND-1
C     IF(((L-1)/LPHNU)*LPHNU.NE.(L-1) .AND. L.NE.ND) GOTO 1
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1
      PRINT 13, (HNU(LPT(I)),I=1,NPRPT)
   13 FORMAT (8(1PE10.3))
C      PRINT 13,L,4.*HNU(L),
C     &          (4.*HNU(L)*RADIUS(L)*RADIUS(L)*PISIG)**0.25
C  13 FORMAT (20X,I10,1PE20.5,0PF15.0)
C   1 CONTINUE
      RETURN
      END subroutine
      end module