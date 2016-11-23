      module mod_inv

      contains

      subroutine inv_lapack(n, a)

      integer, intent(in) :: n

      real*8, intent(inout), dimension(n, n) :: a

      real*8,  dimension(n) :: work ! work array for LAPACK
      integer, dimension(n) :: ipiv   ! pivot indices
      integer               :: info

!     External procedures defined in LAPACK
      external DGETRF
      external DGETRI

!     Store A in Ainv to prevent it from being overwritten by LAPACK
!      Ainv = A

!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.
      call DGETRF(n, n, a, n, ipiv, info)

      if (info .ne. 0) stop 'Matrix is numerically singular!'

!     DGETRI computes the inverse of a matrix using the LU factorization
!     computed by DGETRF.
      call DGETRI(n, a, n, ipiv, work, n, info)

      if (info .ne. 0) stop 'Matrix inversion failed!'

      return

      end subroutine inv_lapack

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

C        SUCHE MAXIMALES MATRIXELEMENT
         AMAX = 0.
         ABSAMAX = 0.

         do 30 J = K, N

            L = N + 1 - K

            ISAMAX = 1

            do LL = 1, L

               IF (abs(A(K+LL-1,j)).gt.abs(A(K+isamax-1,j))) isamax = LL

            enddo

            imax = isamax + K - 1

            if (ABSAMAX .gt. abs(A(IMAX, J))) goto 30

            AMAX =    A(IMAX, J)
            ABSAMAX = ABS(AMAX)

            IK(K) = IMAX
            JK(K) = J

   30    continue
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

      END SUBROUTINE inv

      end module mod_inv
