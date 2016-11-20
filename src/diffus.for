      module MOD_DIFFUS
      contains
      SUBROUTINE DIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR)
      use MOD_BNUE
C***  GIVES THE PLANCK FUNCTION, BCORE, AND ITS RADIUS-DERIVATIVE, DBDR, AT
C***  THE INNER BOUNDARY FROM THE GIVEN TEMPERATURE STRATIFICATION.
C***  IN DIFFUSION APPROXIMATION, THEN THE INCIDENT INTENSITY WILL BE
C***      IPLUS = BCORE + DBDR * Z / X
C***  Z = MUE, X = OPACITY
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 ,intent(out):: BCORE,DBDR
      integer,intent(in) :: ND
      real*8 ,intent(in) :: XLAM,T(ND),RADIUS(ND)
      PARAMETER ( ONE = 1.D+0 )
      BCORE = BNUE(XLAM,T(ND))
      BNDM  = BNUE(XLAM,T(ND-1))
      DBDR=(BCORE-BNDM)/(RADIUS(ND-1)-one)
      RETURN
      END subroutine
      end module