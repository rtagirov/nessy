      module MOD_INVTRI
      contains
C**********  MODULNAME: INVTRI    ******* 24/03/87  21.15.14.******    41 KARTEN
      PURE SUBROUTINE INVTRI (A,B,C,Q,N)
C TRIDIAGONALE MATRIX   -A(L)   B(L)  -C(L)
C  RECHTE SEITE Q(L)
C  LOESUNG AUF VEKTOR Q(L)
C  ACHTUNG -- AUCH C(L) WIRD VERAENDERT --
!      IMPLICIT REAL*8(A-H,O-Z)
	implicit none
	!global variables
      integer,intent(in) :: N
      real*8,dimension(N),intent(in) :: A,B
	real*8,intent(inout),dimension(N) :: C,Q
	!local variables
	integer :: I,L,NM
	real*8  :: AI,CI,HI,QI, H
!      DIMENSION A(N),B(N),C(N),Q(N)
c      DIMENSION AHELP(200),BHELP(200)
c      IF (N .GT. 200) then
c         write (6,*) ' INVTRI dimension insufficient'
c         STOP 'ERROR'
c      endif
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
C**  THE BACKWARD ELIMINATION MAY BE SPEEDED UP BY CRAY VECTOR ROUTINE FOLR
c      DO 4 L=2,N
c      AHELP(L)=-C(N+1-L)
c    4 BHELP(L)= Q(N+1-L)
c      BHELP(1)=Q(N)
c      CALL FOLR(N,AHELP,1,BHELP,1)
c      DO 5 L=1,N
c    5 Q(N+1-L)=BHELP(L)
c      RETURN
C***  DEAD BRANCH: RECURSION BY HAND (CAN NOT BE AUTO-VECTORIZED)
c     activated for PC version
      DO 2 I=1,NM
      L=N-I
      QI=Q(L)+C(L)*QI
    2 Q(L)=QI
      RETURN
      END  subroutine
      end module
