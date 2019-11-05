      module MOD_INVTRI

      contains

      PURE SUBROUTINE INVTRI(A, B, C, Q, N)
!      SUBROUTINE INVTRI(A, B, C, Q, N)

!     TRIDIAGONALE MATRIX   -A(L)   B(L)  -C(L)
!     RECHTE SEITE Q(L)
!     LOESUNG AUF VEKTOR Q(L)
!     ACHTUNG -- AUCH C(L) WIRD VERAENDERT --

      implicit none

!     global variables
      integer, intent(in) ::                 N
      real*8, dimension(N), intent(in) ::    A, B
      real*8, intent(inout), dimension(N) :: C, Q

!     local variables
      integer :: I, L, NM
      real*8  :: AI, CI, HI, QI, H

      CI = C(1) / B(1)
      C(1) = CI
      QI = Q(1) / B(1)
      Q(1) = QI
      NM = N - 1

      if (N .eq. 1) return
      if (N .eq. 2) goto 3

      do I = 2, NM

         AI = A(I)
         H = B(I) - AI * CI
         CI = C(I) / H
         C(I) = CI
         QI = (Q(I) + QI * AI) / H
         Q(I) = QI

!         write(*, '(a,3(1x,e15.7))'), 'invtri 1: ', B(I), AI, H

      enddo

!      stop

    3 QI = (Q(N) + QI * A(N)) / (B(N) - CI * A(N))
      Q(N) = QI

      do I = 1, NM

         L = N - I
         QI = Q(L) + C(L) * QI
         Q(L) = QI

!         write(*, *), 'invtri 2:', Q(L)

      enddo

!      stop

      return

      end subroutine

      end module
