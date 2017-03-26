      MODULE MOD_FLGRID

      CONTAINS

      SUBROUTINE FLGRID(NFL, PHI, PWEIGHT, DELTAX)

!     FREQUENCY SCALE ETC. FOR THE CMF LINE CALCULATIONS
!     NFL          NUMBER OF FREQUENCY POINTS (ORIGINAL POINTS)
!     XMAX         MAXIMUM FREQUENCY
!     XK           FREQUENCY POINTS (FALLING, DOPPLER UNITS) -  NOT STORED
!     PHI(K)       NORMALIZED LINE PROFILE
!     PWEIGHT(K)   RENORMALIZED INTEGRATION WEIGHTS

      IMPLICIT REAL*8(A - H, O - Z)

      real*8, allocatable, dimension(:), intent(out) :: phi, pweight

!     WPI = SQRT(PI)
      DATA WPI /1.7724538509055D0/
     
      XMAX = 6.7D0
      NFL = 61

      if (allocated(phi))     deallocate(phi)
      if (allocated(pweight)) deallocate(pweight)

      allocate(phi(nfl), pweight(nfl))
     
      IF (NFL .GT. NFLDIM) STOP "FLGRID: DIMENSION IS TOO SMALL"

      DELTAX = 2.D0 * XMAX / FLOAT(NFL - 1)
      WS = 0.0D0

      DO K = 1, NFL

         XK = XMAX - (K - 1) * DELTAX
 
         PWEIGHT(K) = EXP(-XK * XK)

         WS = WS + PWEIGHT(K)

         PHI(K) = EXP(-XK * XK) / WPI

      ENDDO
     
!     RENORMALIZATION
      DO K = 1, NFL; PWEIGHT(K) = PWEIGHT(K) / WS; ENDDO

      RETURN

      END SUBROUTINE

      END MODULE
