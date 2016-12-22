      module MATOPER

c     the module contains matrix operations necessary for broyden algebra

      contains

      SUBROUTINE VMF(V2, V1, A, N)

!     MULTIPLICATION VECTOR = VECTOR * MATRIX (FULL) - V2 = V1 * A

      implicit real*8(a - h, o - z)

      DIMENSION V1(N), V2(N), A(N, N)

      DO 1 J = 1, N

         SUM = 0.0d0

         DO 2 I = 1, N

    2    SUM = SUM + V1(I) * A(I, J)

    1 V2(J) = SUM

      return

      end subroutine

      SUBROUTINE VSUB(A, B, N)

      implicit real*8(a - h, o - z)

      DIMENSION A(N), B(N)

      do I = 1, N; A(I) = A(I) - B(I); enddo

      RETURN

      END subroutine

      SUBROUTINE VMV(VSUM, V1, V2, N)
!     MULTIPLICATION RESULT = VECTOR * VEKTOR  --  SUM = V1 * V2^T

      implicit real*8(a - h, o - z)

      DIMENSION V1(N), V2(N)

      VSUM = 0.0d0

      do i = 1, n; VSUM = VSUM + V1(I) * V2(I); enddo

      return

      end subroutine

      SUBROUTINE VMT(V2, A, V1, N)
!     MULTIPLICATION VECTOR = MATRIX (FULL) * VEKTOR TRANSFORMIERT --  V2 = A * V1^T

      implicit real*8(a - h, o - z)

      DIMENSION V1(N), V2(N), A(N, N)

      DO 1 J = 1, N

         VSUM = .0

         DO 2 I = 1, N

            VSUM = VSUM + V1(I) * A(J, I)

    2    CONTINUE

         V2(J) = VSUM

    1 CONTINUE

      RETURN

      END subroutine

      SUBROUTINE VMD(A, V1, V2, N)
!     MULTIPLICATION MATRIX = VECTOR * VEKTOR, DYADIC PRODUCT -- A = V1^T * V2

      implicit real*8(a - h, o - z)

      DIMENSION V1(N), V2(N), A(N, N)

      DO 1 J = 1 ,N

         DO 2 I = 1, N

            A(J, I) = V1(J) * V2(I)

    2    CONTINUE

    1 CONTINUE

      RETURN

      END subroutine

      SUBROUTINE VDIFF(A1, VDI, N)
!     DIVISION MATRIX = MATRIX / FLOAT  --  A1 = A2 / VDIFF

      implicit real*8(a - h, o - z)

      DIMENSION A1(N, N)

      IF (VDI .EQ. .0) GOTO 3

      DO 2 J = 1, N

         DO 1 I = 1, N

         A1(J, I) = A1(J, I) / VDI

    1    CONTINUE

    2 CONTINUE

    3 CONTINUE

      RETURN

      END subroutine

      SUBROUTINE VADD (A, B, N)
!     VECTOR ADDITION  A = A + B

      implicit real*8(a - h, o - z)

      DIMENSION A(N), B(N)

      do 1 I = 1, N

    1 A(I) = A(I) + B(I)

      return

      end subroutine

      SUBROUTINE VADDM (A1, A2, N)
!     ADDITION OF MATRIX A1 = A1 + A2

      implicit real*8(a - h, o - z)

      DIMENSION A1(N, N), A2(N, N)

      DO 1 J = 1, N

         DO 2 I = 1, N

            A1(J, I) = A1(J, I) + A2(J, I)

    2    CONTINUE

    1 CONTINUE

      return

      end subroutine

      SUBROUTINE ACOPY (A1, A2, N)

      implicit real*8(a - h, o - z)

      dimension A1(N, N), A2(N, N)

      do J = 1, N

         do I = 1, N

            A1(J, I) = A2(J, I)

         enddo

      enddo

      return

      end subroutine

      subroutine inv_lapack(n, a)

      integer, intent(in) :: n

      real*8, intent(inout), dimension(n, n) :: a

      real*8,  dimension(n) :: work ! work array for LAPACK
      integer, dimension(n) :: ipiv ! pivot indices
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

      RETURN

      END SUBROUTINE inv

      SUBROUTINE MVMD(BX, B, C, JMAX, JMM)
C***  MATRIX (VOLL)  B  *  MATRIX (DIAGONAL)  C
C***  ERGEBNIS-MATRIX  BX(VOLL)
C***  AKTUELLES FORMAT BX(JMAX,JMM)=B(JMAX,JMAX)*C(JMAX,JMM)
C***  WOBEI DIE UEBERZAEHLIGEN ZEILEN DER DIAGONALMATRIX C VERSCHWINDEN
      implicit real*8(a - h, o - z)

      DIMENSION BX(JMAX, JMAX), B(JMAX, JMAX), C(JMAX)
      DO 1 K = 1, JMM
      CK = C(K)
      DO 1 I = 1, JMAX
    1 BX(I, K) = B(I, K) * CK

      RETURN

      end subroutine

      SUBROUTINE MVV(WX, B, W, JMAX, JMM)
!     MATRIX (VOLL)  B  *  VEKTOR W
!     ERGEBNIS-VEKTOR  WX
!     AKTUELLES FORMAT  WX(JMAX) = B(JMAX,JMM) * W(JMM)
      implicit real*8(a - h, o - z)

      DIMENSION WX(JMAX), B(JMAX, JMAX), W(JMAX)
      DO 1 I=1,JMAX
      WXI = .0
      DO 2 K=1,JMM
    2 WXI=WXI+B(I,K)*W(K)
    1 WX(I)=WXI

      RETURN

      end subroutine

      SUBROUTINE MDV(A, W, N)
!     MATRIX A (DIAGONAL)  *  VEKTOR W
!     ERGEBNIS-VEKTOR UEBERSCHREIBT  W

      implicit real*8(a - h, o - z)

      DIMENSION A(N), W(N)
      DO 1 I = 1, N
    1 W(I) = A(I) * W(I)

      RETURN

      end subroutine

      SUBROUTINE MSUB(A, B, N)
!     A = A - B
      implicit real*8(a - h, o - z)

      DIMENSION A(N, N), B(N, N)

      DO 1 I = 1, N
      DO 1 K = 1, N
    1 A(I, K) = A(I, K) - B(I, K)

      RETURN

      end subroutine

      SUBROUTINE MDMV(A, B, N)
!     MATRIX (DIAGONAL)  A  *  MATRIX (VOLL)  B
!     ERGEBNIS-MATRIX UEBERSCHREIBT B
      implicit real*8(a - h, o - z)

      DIMENSION A(N), B(N, N)
      DO 1 I = 1, N
      AI = A(I)
      DO 1 K = 1, N
    1 B(I, K) = B(I, K) * AI

      RETURN

      end subroutine

      SUBROUTINE VMALV(VA, VB, V, Q, LMAX)
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
      implicit real*8(a - h, o - z)

      DIMENSION VA(LMAX), VB(LMAX), V(LMAX), Q(LMAX)
      LZ=LMAX-1
      Q(1)=VB(1)*V(1)

      DO 1 L=2,LZ
      Q(L)=VA(L)*V(L-1)+VB(L)*V(L)
    1 CONTINUE

      Q(LMAX)=VA(LMAX)*V(LZ)

      RETURN

      END subroutine

      end module
