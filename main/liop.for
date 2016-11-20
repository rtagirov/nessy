      MODULE MOD_LIOP

      CONTAINS

      SUBROUTINE LIOP(EINST, WEIGHTI, WEIGHTJ, I, J,
     $                ND, XLAM, ENTOT, POPNUM, RSTAR,
     $                OPAL, ETAL, VDOP,
     $                POPNUM_SECOND_DIM_SIZE)
     
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

      INTEGER, INTENT(IN) :: POPNUM_SECOND_DIM_SIZE

      DIMENSION ENTOT(ND), OPAL(ND), ETAL(ND)

      REAL*8, DIMENSION(ND, POPNUM_SECOND_DIM_SIZE) :: POPNUM

      REAL*8, INTENT(IN) :: ENTOT

!     C3 = 4 * PI / H / C (CGS UNITS)
      DATA C3 /6.3268D16/

!     PI8 = 8 * PI
      DATA PI8 /25.1327D0/

      XLAMCM = XLAM / 1.0D+8

!     DND = DELTA-NUE-DOPPLER (HERTZ)
!     ************************************************************************************
      DND = VDOP / XLAM * 1.0D+13

!     RINAT TAGIROV:
!     VDOP = sqrt(V_therm^2 + V_turb^2), where V_therm is the
!     thermal velocity of atoms and V_turb is the turbulent velocity
!     of separate parcels of gas. This can be inferred from the FLGRID.for procedure,
!     formula for the Doppler line profile in page 205 of Lang's "Astrophysical formulas"
!     and the way VDOP enters this routine. While it is very simple to calculate
!     VDOP at any height having the model of the atmosphere in HMINUS
!     it is a free parameter set it the CARDS file.
!     It has to be free because of some computational obstacles arising in case
!     when VDOP is depth dependent (for details refer to Werner cuz I never quite figured
!     out the reason properly). However, one of the possible reasons could be that 
!     in ETL.FOR the velocity field and its gradient are divided by VDOP
!     for the purposes of their conversion into dimensionless units.
!     Velocity field has to be monotonic and dividing it by a function of height
!     could obviate the monotonicity. Therefore VDOP has to be constant, hence it has to
!     be a free parameter. Why one needs to introduce the dimensionless velocity field
!     is another question which I don't have an answer to.
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

      END MODULE
