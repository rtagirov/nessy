      module storextr

      implicit none

      contains

      subroutine extrxjc(XJCREA, XJC, EDDREA, EDDI, nd, nf, K)

      integer, intent(in) :: nd, nf, K

      real*8,  intent(in)  :: xjcrea(nd, nf), eddrea(3, nd, nf)

      real*8,  intent(out) :: xjc(nd), eddi(3, nd)

      xjc(1 : ND) =     xjcrea(1 : ND, K)

      eddi(1, 1 : ND) = eddrea(1, 1 : ND, K)
      eddi(2, 1 : ND) = eddrea(2, 1 : ND, K)
      eddi(3, 1 : ND) = eddrea(3, 1 : ND, K)

      return

      end subroutine

      subroutine storxjc(XJCREA, XJC, EDDREA, EDDI, nd, nf, K)

      integer, intent(in) :: nd, nf, K

      real*8,  intent(in) :: xjc(nd), eddi(3, nd)

      real*8,  intent(out)  :: xjcrea(nd, nf), eddrea(3, nd, nf)

      xjcrea(1 : ND, K) =    xjc(1 : ND)

      eddrea(1, 1 : ND, K) = eddi(1, 1 : ND)
      eddrea(2, 1 : ND, K) = eddi(2, 1 : ND)
      eddrea(3, 1 : ND, K) = eddi(3, 1 : ND)

      return

      end subroutine

      subroutine storxjl(XJL, XJLMEAN, ND, LASTIND, IND, LLO, LO)

      integer, intent(in) ::            ND, IND, LASTIND

      real*8, dimension(ND, LASTIND) :: XJL

      real*8, dimension(ND) ::          XJLMEAN

!     the Local approximate lambda-Operator for the line with index IND
      real*8, dimension(ND) ::          LO

!     the Local approximate lambda-Operator for all lines (Overall)
      real*8, dimension(ND, LASTIND) :: LLO

      XJL(1 : ND, IND) = XJLMEAN(1 : ND)

      LLO(1 : ND, IND) = LO(1 : ND)

      return

      end subroutine

      end module
