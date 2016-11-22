      module mod_inv

      contains

      subroutine inv(n, a)

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

      end subroutine inv

      end module mod_inv
