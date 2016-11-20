      MODULE MOD_HYD_ION_STEP
      CONTAINS


      SUBROUTINE HYD_ION_STEP(DepthIndex, IonizationDepthIndex, 
     $ IonPopul, ElecConc)


      INTEGER, INTENT(IN) :: DepthIndex, IonizationDepthIndex

      REAL*8, DIMENSION(12), INTENT(INOUT) :: IonPopul
      REAL*8,                INTENT(INOUT) :: ElecConc

      REAL*8  :: TotNeutHydConc, TotHydConc
      REAL*8  :: HydMinusConc,   ProtConc


      TotNeutHydConc = SUM(IonPopul(1:11))
      TotHydConc =     SUM(IonPopul(1:12))

      HydMinusConc =   IonPopul(1)
      ProtConc =       IonPopul(12)

      IF (DepthIndex .LE. IonizationDepthIndex) THEN

         IonPopul(1:11) = 0.
         IonPopul(12) =   TotHydConc

      ELSE

         DO i = 1, 11

            IonPopul(i) = IonPopul(i) * TotHydConc / TotNeutHydConc

         ENDDO

         IonPopul(12) = 0.

      ENDIF

      ElecConc = ElecConc + (IonPopul(12) - ProtConc) -
     $                      (IonPopul(1)  - HydMinusConc)

      RETURN

      END SUBROUTINE


      END MODULE
