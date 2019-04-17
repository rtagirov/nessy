      MODULE MOD_LIOP

      CONTAINS

      SUBROUTINE LIOP_RTE(EINST, WEIGHTI, WEIGHTJ, I, J,
     $                    ND, XLAM, ENTOT, POPNUM, RSTAR,
     $                    OPAL, ETAL, VDOP, N)
     
!     CALCULATES THE LINE OPACITY OPAL AND THE EMISSIVITY ETAL FOR ONE LINE
!     AT ALL DEPTH POINTS
!     PHYSICAL DIMENSION OF OPAL : PER RSTAR AND PER DELTA-NUE-DOPPLER
!     (THE LATTER WILL BE CANCELLED OUT BY THE NORMALISED PROFILE FUNCTION)
!     GIVEN QUANTITIES : EINST = EINSTEIN COEFFICIENT A (UP - LOW, PER SECOND)
!                        WEIGHTI, WEIGHTJ = STATISTICAL WEIGHTS (UP, LOW)
!                        POPNUM(L, J) = RELATIVE POPULATION NUMBERS
!                        ENTOT(L) = TOTAL NUMBER DENSITY

      IMPLICIT REAL*8(A - H, O - Z)

!     When liop.for is called from etl.for ND = ND, when called from other routines ND = 1
      INTEGER, INTENT(IN) :: ND

      INTEGER, INTENT(IN) :: N

      DIMENSION ENTOT(ND), OPAL(ND), ETAL(ND)

      REAL*8, DIMENSION(ND, N) :: POPNUM

      REAL*8, INTENT(IN) :: ENTOT

!     C3 = 4 * PI / H / C (CGS UNITS)
      DATA C3 /6.3268D16/

!     PI8 = 8 * PI
      DATA PI8 /25.1327D0/

      XLAMCM = XLAM / 1.0D+8

!     DND = DELTA-NUE-DOPPLER (HERTZ)
!     ************************************************************************************
      DND = VDOP / XLAM * 1.0D+13

      print*, 'liop_rte VDOP: ', VDOP

!     RINAT TAGIROV:
!     VDOP = sqrt(V_therm^2 + V_turb^2), where V_therm is the
!     thermal velocity of atoms and V_turb is the turbulent velocity
!     of separate parcels of gas.
!     It is a free parameter set in the CARDS file.
!     One of the reasons for it to be free is that
!     in ETL.FOR the velocity field and its gradient are divided by VDOP
!     for them to be converted into dimensionless units.
!     Velocity field has to be monotonic and dividing it by a function of height
!     could violate the monotonicity.
!     ************************************************************************************

      EMINDU = EINST * XLAMCM * XLAMCM / PI8 * RSTAR
      ABSORP = EMINDU * WEIGHTJ / WEIGHTI
      EMSPON = EINST / XLAMCM / C3 * RSTAR

      DO L = 1, ND

         ENI = POPNUM(L, I) * ENTOT(L)
         ENJ = POPNUM(L, J) * ENTOT(L)

         OPAL(L) = (ENI * ABSORP - ENJ * EMINDU) / DND

         IF (OPAL(L) .LT. 0.0D0) OPAL(L) = 0.0D0

         ETAL(L) = ENJ * EMSPON / DND

      ENDDO

      RETURN

      END SUBROUTINE


      SUBROUTINE LIOP_SBE(EINST, WEIGHTI, WEIGHTJ, I, J,
     $                    XLAM, ENTOTL, POPNUM, RSTAR,
     $                    OPAL, ETAL, VDOP, N)
     
!     CALCULATES THE LINE OPACITY OPAL AND THE EMISSIVITY ETAL FOR ONE LINE
!     AT ALL DEPTH POINTS
!     PHYSICAL DIMENSION OF OPAL : PER RSTAR AND PER DELTA-NUE-DOPPLER
!     (THE LATTER WILL BE CANCELLED OUT BY THE NORMALISED PROFILE FUNCTION)
!     GIVEN QUANTITIES : EINST = EINSTEIN COEFFICIENT A (UP - LOW, PER SECOND)
!                        WEIGHTI, WEIGHTJ = STATISTICAL WEIGHTS (UP, LOW)
!                        POPNUM(L, J) = RELATIVE POPULATION NUMBERS
!                        ENTOT(L) = TOTAL NUMBER DENSITY

      IMPLICIT REAL*8(A - H, O - Z)

!     When liop.for is called from etl.for ND = ND, when called from other routines ND = 1
      INTEGER, INTENT(IN) :: N

      REAL*8, DIMENSION(N) :: POPNUM

      REAL*8, INTENT(IN) :: ENTOTL

      real*8 :: opal, etal

!     C3 = 4 * PI / H / C (CGS UNITS)
      DATA C3 /6.3268D16/

!     PI8 = 8 * PI
      DATA PI8 /25.1327D0/

      XLAMCM = XLAM / 1.0D+8

!     DND = DELTA-NUE-DOPPLER (HERTZ)
!     ************************************************************************************
      DND = VDOP / XLAM * 1.0D+13

      print*, 'liop_sbe VDOP: ', VDOP

!     RINAT TAGIROV:
!     VDOP = sqrt(V_therm^2 + V_turb^2), where V_therm is the
!     thermal velocity of atoms and V_turb is the turbulent velocity
!     of separate parcels of gas.
!     It is a free parameter set in the CARDS file.
!     One of the reasons for it to be free is that
!     in ETL.FOR the velocity field and its gradient are divided by VDOP
!     for them to be converted into dimensionless units.
!     Velocity field has to be monotonic and dividing it by a function of height
!     could violate the monotonicity.
!     ************************************************************************************

      EMINDU = EINST * XLAMCM * XLAMCM / PI8 * RSTAR
      ABSORP = EMINDU * WEIGHTJ / WEIGHTI
      EMSPON = EINST / XLAMCM / C3 * RSTAR

      ENI = POPNUM(I) * ENTOTL
      ENJ = POPNUM(J) * ENTOTL

      OPAL = (ENI * ABSORP - ENJ * EMINDU) / DND

      ETAL = ENJ * EMSPON / DND

      IF (OPAL .LT. 0.0D0) OPAL = 0.0D0

      RETURN

      END SUBROUTINE

      END MODULE
