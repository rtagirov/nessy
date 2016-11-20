      MODULE FIOSS_AUX

      CONTAINS

      SUBROUTINE READ_NLTETRAPOP(TI, TPL, TPU, TDL, TDU)

      USE COMMON_BLOCK
      USE FILE_OPERATIONS

      IMPLICIT NONE

      INTEGER, DIMENSION(NTC), INTENT(IN) ::               TI

      REAL*8, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: TPL, TPU, TDL, TDU

      CHARACTER (LEN = 32) :: FUDGE

      INTEGER :: I, J, L, DEPTH_POINTS_NUM

      REAL*8 :: PL, PU, DL, DU

      DEPTH_POINTS_NUM = NUM_OF_LINES(ATM_MOD_FILE)

      ALLOCATE(TPL(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TPU(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TDL(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TDU(DEPTH_POINTS_NUM, NTC))

      TPL(:, :) = 0.0D0; TPU(:, :) = 0.0D0; TDL(:, :) = 0.0D0; TDU(:, :) = 0.0D0
 
      DO I = 1, NTC

         OPEN(773, FILE = NTP_FILE)

         IF (TI(I) .NE. 1) THEN

             DO J = 1, (TI(I) - 1) * DEPTH_POINTS_NUM; READ(773, '(A66)') FUDGE; ENDDO

         ENDIF

         DO L = 1, DEPTH_POINTS_NUM

            READ(773, '(E15.7,3(2x,E15.7))') PL, PU, DL, DU

            TPL(L, I) = PL; TPU(L, I) = PU; TDL(L, I) = DL; TDU(L, I) = DU

            WRITE(*, '(A,4(1x,E15.7))') 'TEST READ NLTE:', TPL(L, I), TPU(L, I), TDL(L, I), TDU(L, I)

         ENDDO

         CLOSE(773)

      ENDDO

      END SUBROUTINE


      SUBROUTINE READ_NLTE_TRA(NFE, NUMTR, EID, CH, WL, WU, AUPLOW, WAVTR)

      USE FILE_OPERATIONS

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) ::                            NFE

      INTEGER, INTENT(OUT) ::                            NUMTR

      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: EID, WL, WU, CH

      REAL*8, ALLOCATABLE, DIMENSION(:), INTENT(OUT) ::  WAVTR, AUPLOW

      INTEGER ::                                         I, EI, GU, GL, CHARG

      CHARACTER (LEN = 10) ::                            LL, LU

      REAL*8 ::                                          WAVT, AUL

      INQUIRE(FILE = NTW_FILE, EXIST = NFE)

      IF (NFE .EQ. .FALSE.) RETURN

      NUMTR = NUM_OF_LINES(NTW_FILE)

      ALLOCATE(AUPLOW(NUMTR))
      ALLOCATE(WAVTR(NUMTR))
      ALLOCATE(EID(NUMTR))    ! ELEMENT IDENTIFICATION, I.E. ITS NUMBER IN THE PERIODICAL TABLE
      ALLOCATE(CH(NUMTR))
      ALLOCATE(WL(NUMTR))
      ALLOCATE(WU(NUMTR))

      OPEN(132, FILE = NTW_FILE)

      DO I = 1, NUMTR

         READ(132, '(2(A10,2x),I2,2x,I2,2(2x,I3),2(2x,E15.7))') LL, LU, EI, CHARG, GL, GU, AUL, WAVT

         WAVTR(I) = WAVT; EID(I) = EI; CH(I) = CHARG; WL(I) = GL; WU(I) = GU; AUPLOW(I) = AUL

         WRITE(*, '(A,I2,3(2x,I3),2(2x,E15.7))') 'NLTE TRA:', EID(I), CH(I), WL(I), WU(I), AUPLOW(I), WAVTR(I)

      ENDDO

      END SUBROUTINE


      FUNCTION airlambda(vaclambda)
!    translate vacuum to the airlambda lambda in A                             
                                                                                
      IMPLICIT NONE
      real*8 vaclambda, airlambda, sig, n

      sig=1.d4/(vaclambda)
      n=1.d0+6.4328d-5+(2.94981d-2)/(146.-sig**2.)+(2.554d-4)/
     $ (41-sig**2.)
      airlambda=vaclambda/n
      return
      END FUNCTION


      END MODULE FIOSS_AUX
