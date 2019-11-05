      module storextr

      implicit none

      contains

      subroutine extrxjc(xjc2, xjc1, eddi3, eddi2, nd, nf, k)

      integer, intent(in) :: nd, nf, k

      real*8,  intent(in)  :: xjc2(nd, nf), eddi3(3, nd, nf)

      real*8,  intent(out) :: xjc1(nd), eddi2(3, nd)

      xjc1(1 : ND) =     xjc2(1 : ND, k)

      eddi2(1, 1 : ND) = eddi3(1, 1 : ND, k)
      eddi2(2, 1 : ND) = eddi3(2, 1 : ND, k)
      eddi2(3, 1 : ND) = eddi3(3, 1 : ND, k)

      return

      end subroutine

      subroutine storxjc(xjc2, xjc1, eddi3, eddi2, nd, nf, k)

      integer, intent(in) :: nd, nf, k

      real*8,  intent(in) :: xjc1(nd), eddi2(3, nd)

      real*8,  intent(out)  :: xjc2(nd, nf), eddi3(3, nd, nf)

      xjc2(1 : ND, k) =    xjc1(1 : ND)

      eddi3(1, 1 : ND, k) = eddi2(1, 1 : ND)
      eddi3(2, 1 : ND, k) = eddi2(2, 1 : ND)
      eddi3(3, 1 : ND, k) = eddi2(3, 1 : ND)

      return

      end subroutine

      subroutine storxjl(xjl2, xjl1, nd, lastind, ind, lo2, lo1)

      integer, intent(in) ::            nd, ind, lastind

      real*8, dimension(nd, lastind) :: xjl2

      real*8, dimension(nd) ::          xjl1

!     the Local approximate lambda-Operator for the line with index IND
      real*8, dimension(nd) ::          lo1

!     the Local approximate lambda-Operator for all lines (Overall)
      real*8, dimension(nd, lastind) :: lo2

      xjl2(1 : ND, IND) = xjl1(1 : ND)

      lo2(1 : ND, IND) =  lo1(1 : ND)

      return

      end subroutine

      end module
