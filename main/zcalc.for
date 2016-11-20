      module MOD_ZCALC
      contains
      SUBROUTINE ZCALC (RADIUS,P,Z,ND,NP)
      implicit none
C***  THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN Z
      real*8,intent(in)  :: RADIUS(ND),P(NP)
      real*8,intent(out) :: Z(ND*NP)
      integer,intent(in) :: ND,NP
      integer            :: I,J,L,JMAX
      real*8             :: RR
      DO L=1,ND
        RR=RADIUS(L)*RADIUS(L)
        JMAX=NP+1-L
        DO J=1,JMAX
          I=(J-1)*ND+L
          Z(I)=SQRT(RR-P(J)**2)
        ENDDO
      ENDDO
      RETURN
      END subroutine
      end module