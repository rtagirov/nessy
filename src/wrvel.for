      MODULE MOD_WRVEL
      CONTAINS

      FUNCTION WRVEL(R)

      IMPLICIT NONE

      real*8, intent(in) :: R
      real*8 ::             WRVEL
      real*8 ::             RCON,VPAR1,VPAR2,BETA,VFINAL,VMIN,HSCALE

      COMMON/VELPAR/        VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE

      IF (R .LE. 1.d0) GOTO 2
      IF (R .LT. RCON) GOTO 3

      WRVEL = VPAR1 * (1.d0 - VPAR2 / R)**BETA

      IF (WRVEL .GT. VFINAL) WRVEL = VFINAL

      RETURN

    2 WRVEL = VMIN

      RETURN

    3 WRVEL = VMIN / R / R * EXP((1.d0 - 1.d0 / R) / HSCALE)

      RETURN

      END FUNCTION

      END MODULE
