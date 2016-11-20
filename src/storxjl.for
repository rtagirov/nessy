      MODULE MOD_STORXJL

      CONTAINS

      SUBROUTINE STORXJL(XJL, XJLMEAN, ND, LASTIND, IND, LOO, LO)

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::            ND, IND, LASTIND

      REAL*8, DIMENSION(ND, LASTIND) :: XJL

      REAL*8, DIMENSION(ND) ::          XJLMEAN

!     the Local approximate lambda-Operator for the line with index IND
      REAL*8, DIMENSION(ND) ::          LO

!     the Local approximate lambda-Operator for all lines (Overall)
      REAL*8, DIMENSION(ND, LASTIND) :: LOO

      XJL(1 : ND, IND) = XJLMEAN(1 : ND)

      LOO(1 : ND, IND) = LO(1 : ND)

      RETURN

      END SUBROUTINE

      END MODULE
