      module MOD_INVTRI

      contains

      PURE SUBROUTINE INVTRI(A, B, C, Q, N)

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

      CI=C(1)/B(1)
      C(1)=CI
      QI=Q(1)/B(1)
      Q(1)=QI
      NM=N-1
      IF (N.EQ.1) RETURN
      IF (N.EQ.2) GOTO 3
      DO 1 I=2,NM
      AI=A(I)
      H=B(I)-AI*CI
      CI=C(I)/H
      C(I)=CI
      QI=(Q(I)+QI*AI)/H
    1 Q(I)=QI
    3 QI=(Q(N)+QI*A(N))/(B(N)-CI*A(N))
      Q(N)=QI

      DO 2 I=1,NM
      L=N-I
      QI=Q(L)+C(L)*QI
    2 Q(L)=QI

      RETURN

      END  subroutine

      end module
