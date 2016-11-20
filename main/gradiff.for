      module MOD_GRADIFF
      contains
C**********  MODULNAME: GRADIFF   ******* 24/03/87  21.12.18.******    13 KARTEN
      SUBROUTINE GRADIFF  (ND,VELO,GRADI,RADIUS)
C***  FOR THE VELOCITY FIELD GIVEN BY VECTOR VELO(L), THE GRADIENTS GRADI(L)
C***  ARE COMPUTED BY LINEAR INTERPOLATION BETWEEN THE NEIGHBORING POINTS

      IMPLICIT REAL*8(A-H,O-Z)
      real*8,intent(out) ::GRADI
      real*8,intent(in) ::VELO,RADIUS
      integer,intent(in)::ND
      DIMENSION VELO(ND),GRADI(ND),RADIUS(ND)
      NDM=ND-1
      forall(L=2:NDM)
        GRADI(L)=(VELO(L+1)-VELO(L-1))/(RADIUS(L+1)-RADIUS(L-1))
      endforall
      GRADI(1)=GRADI(2)
      GRADI(ND)=GRADI(NDM)

      RETURN
      END subroutine
      end module