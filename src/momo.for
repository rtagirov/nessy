      module MOD_MOMO

      contains

      SUBROUTINE MOMO (OPA,ETA,THOMSON,EDDI,R,XJC,ND)

!     SOLUTION OF THE MOMENT EQUATION
!     Calculates XJC 
!     input: EDDI, OPA, ETA, THOMSO, R
!     INVTRI inverts the triangular matrix [A,B,C], changes C

      use MOD_INVTRI

      implicit none

      integer, intent(in) ::                  ND

      real*8, intent(in), dimension(ND) ::    ETA, OPA, THOMSON, R
      real*8, intent(in), dimension(3, ND) :: EDDI

      real*8, intent(out), dimension(ND) ::   XJC

      real*8, dimension(ND) :: A, B, C, W
	
      real*8 ::  DA, DC, DR, DX, FL, FM, FP, QL, QM, QP, H, HPLUS
      real*8 ::  RL, X
      integer :: L, NDM
	
      real*8, parameter :: one = 1.0d0, two = 2.d0, four = 4.0d0
     
!     OUTER BOUNDARY
      FL=EDDI(1,1)
      QL=EDDI(2,1)
      H =EDDI(3,1)
      FP=EDDI(1,2)
      QP=EDDI(2,2)
      RL=R(1)
      X=OPA(1)
      DX=(QL+QP)*(X+OPA(2))*(RL-R(2))/four
      B(1)=two*QL*QL*FL*X/DX/DX+two*QL*X*H/DX+X*(one-THOMSON(1))
      C(1)=two*QL*QP*FP*X/DX/DX
      W(1)=RL*RL*ETA(1)
     
!     NON-BOUNDARY POINTS
      NDM = ND - 1

      DO L=2,ND-1
         FL=EDDI(1,L)
         QL=EDDI(2,L)
         FP=EDDI(1,L+1)
         QP=EDDI(2,L+1)
         FM=EDDI(1,L-1)
         QM=EDDI(2,L-1)
         X=OPA(L)
         RL=R(L)
         DR=(R(L-1)-R(L+1))/two
         DA=DR*(R(L-1)-RL)*(QM+QL)*(OPA(L-1)+X)/four
         DC=DR*(RL-R(L+1))*(QP+QL)*(OPA(L+1)+X)/four
         A(L)=QM*FM/DA
         C(L)=QP*FP/DC
         B(L)=QL*FL*(one/DA+one/DC)+X*(one-THOMSON(L))
         W(L)=RL*RL*ETA(L)
      ENDDO
     
!     INNER BOUNDARY
      FL=EDDI(1,ND)
      QL=EDDI(2,ND)
      H =EDDI(3,ND)
      HPLUS=EDDI(3,ND-1)
      FM=EDDI(1,ND-1)
      QM=EDDI(2,ND-1)
      RL=R(ND)
      X=OPA(ND)
      DX=(QM+QL)*(X+OPA(ND-1))*(R(ND-1)-R(ND))/four
      A(ND)=two*QL*QM*FM*X/DX/DX
      B(ND)=two*QL*QL*FL*X/DX/DX+two*QL*X*H/DX+X*(one-THOMSON(ND))
      W(ND)=RL*RL*(ETA(ND)+two*QL*X*HPLUS/DX)
     
      call invtri(A, B, C, W, ND)

      XJC(1 : ND) = W(1 : ND) / R(1 : ND) / R(1 : ND)

      return

      end subroutine

      end module
