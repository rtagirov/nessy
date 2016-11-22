      module MOD_INV

      contains

      SUBROUTINE INV(N, A)

C******************************************************************************
C***  MATRIX INVERSION (CRAY FORTRAN)
C***  MATRIX A
C***  N = RANK OF SUBMATRIX (LEFT UPPER BLOCK) TO BE INVERTED
C******************************************************************************

      USE MOD_ERROR

      implicit real*8(a - h, o - z)

      integer, intent(in) :: N

      real*8, intent(inout) :: A(N, N)

      integer :: IK(N), JK(N)

      DO 100 K = 1, N
C
C     SUCHE MAXIMALES MATRIXELEMENT
!      AMAX=0.
!      ABSAMAX=0.

!         do J = K, N
!        L=N+1-K
!        ISAMAX=1
!        DO LL=1,L
c           PRINT *, 'N=',N, 'J=',J,'LL=',LL,'A(LL,J)=',A(LL,J) 
!           IF (abs(A(K+LL-1,j)).gt.abs(A(K+isamax-1,j))) isamax=LL
!        ENDDO
!        imax=isamax+K-1
c*cray      IMAX=ISAMAX(L,A(K,J),1)+K-1
!      IF(ABSAMAX.GT.ABS(A(IMAX,J))) GOTO 30
!      AMAX=A(IMAX,J)
!      ABSAMAX=ABS(AMAX)

            imax = maxloc(abs(A(k, k : n))

            amax = a(imax, j)

            absamax = abs(amax)

            IK(K) = IMAX
            JK(K) = J

!         enddo
C
C***  STOP IN CASE OF SINGULARITY (PIVOT ELEMENT = 0 )
      IF (AMAX .EQ. .0) THEN
        print *,'INV: SINGULARITY DISCOVERED'
        print *,'INV: ',N, AMAX
        CALL ERROR('inv.for: ERROR SINGULARITY DISCOVERED')
      ENDIF
C
C     VERTAUSCHEN DER ZEILEN I UND K
      I=IK(K)
      IF(I.NE.K) THEN
        DO J=1,N
          SAVE =A(K,J)
          A(K,J)=A(I,J)
          A(I,J)=-SAVE
        ENDDO
      ENDIF
C
C     VERTAUSCHEN DER SPALTEN J UND K
      J=JK(K)
      IF(J.NE.K) THEN
        DO I=1,N
          SAVE = A(I,K)
          A(I,K)=A(I,J)
          A(I,J)=-SAVE
        ENDDO
      ENDIF
C
C     DIVISION DER SPALTE K DURCH AMAX
      A(1:N,K)=-A(1:N,K)/AMAX
      A(K,K)=0.
C
C     UMFORMEN DER ZEILEN UND SPALTEN
C***  ELEMENTARE UMFORMUNG:  ZEILE I = ZEILE I + A(I,K) * ZEILE K
      DO 81 I=1,N
      IF (I .EQ. K) GOTO 81
      AIK=A(I,K)
      A(I,1:N)=A(I,1:N)+AIK*A(K,1:N)
   81 CONTINUE
C
C***  SPALTE K: DIVISION DUCH AMAX
      DO 90 J=1,N
      A(K,J)=A(K,J)/AMAX
  90  CONTINUE
C
C***  DIAGONALELEMENT
      A(K,K)=1./AMAX
  100 CONTINUE
C
C
C     ZEILEN UND SPALTEN RUECKTAUSCHOPERATIONEN
      DO 130 L=1,N
      K=N-L+1
      J=IK(K)
      IF(J.LE.K) GOTO 111
C***  VERTAUSCHEN DER SPALTEN J UND K
      DO 110 I=1,N
       SAVE=A(I,K)
      A(I,K)=-A(I,J)
  110 A(I,J)=SAVE
  111 I=JK(K)
      IF(I .LE. K) GOTO 130
C***  VERTAUSCHEN DER ZEILEN I UND K
      DO 120 J=1,N
      SAVE=A(K,J)
      A(K,J)=-A(I,J)
  120 A(I,J)=SAVE
  130 CONTINUE
C
      RETURN
      END SUBROUTINE
      END MODULE
