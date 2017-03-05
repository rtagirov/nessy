      MODULE COMMON_BLOCK

!     HMINUS VARIABLES

      INTEGER ::                                 DPN

      REAL*8 ::                                  CORMAX_ELEC
 
      LOGICAL ::                                 LTE_RUN

      REAL*8, DIMENSION(:, :), ALLOCATABLE ::    XJL_LTE
      REAL*8, DIMENSION(:, :), ALLOCATABLE ::    XJC_LTE
      REAL*8, DIMENSION(:, :), ALLOCATABLE ::    POP_LTE

      REAL*8, DIMENSION(:, :, :), ALLOCATABLE :: ARR, ACR, RBR
      REAL*8, DIMENSION(:, :, :), ALLOCATABLE :: ARR_LTE

      CHARACTER(:), ALLOCATABLE ::               NRRM_FILE, NCRM_FILE, NTRM_FILE

!     the Local approximate lambda-Operator for all lines (Overall)
!******************************************************************
      REAL*8, DIMENSION(:, :), ALLOCATABLE ::    LLO
!******************************************************************

      real*8, dimension(:, :), allocatable ::    tau_line, tau_cont

      character (len = 1), dimension(:, :), allocatable :: damp_line, damp_cont

!      logical ::                                 damp_acc = .true.
      logical ::                                 damp_acc = .false.

      logical ::                                 rayleigh = .true.

      INTEGER ::                                 LAMBDA_ITER

      LOGICAL, DIMENSION(:), ALLOCATABLE ::      NOFILE

      LOGICAL ::                                 CONST_ELEC = .FALSE.

      logical ::                                 each_ali = .false.

      LOGICAL ::                                 VEL_FIELD_FROM_FILE

      REAL*8, DIMENSION(:), ALLOCATABLE ::       HEIGHT

      LOGICAL ::                                 OLDSTART

      INTEGER ::                                 NLINE

!     FIOSS VARIABLES

      INTEGER ::                                 NUMTRA

      INTEGER ::                                 NTC

      INTEGER, ALLOCATABLE, DIMENSION(:) ::      NTI

      INTEGER, ALLOCATABLE, DIMENSION(:) ::      ELID, CHARGE, SWL, SWU

      REAL*8,  ALLOCATABLE, DIMENSION(:) ::      WAVTRA, AUL

      REAL*8,  ALLOCATABLE, DIMENSION(:, :) ::   NTPL, NTPU, NTDL, NTDU

!      REAL*8,  ALLOCATABLE, DIMENSION(:) ::      NLTE_ABS, NLTE_EMI

      END MODULE COMMON_BLOCK
