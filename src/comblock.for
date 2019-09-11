      module common_block

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

      logical, dimension(:, :), allocatable ::   damp_line, damp_cont

      logical ::                                 damp_acc = .true.
!      logical ::                                 damp_acc = .false.

!      logical ::                                 rayleigh = .true.
      logical ::                                 rayleigh = .false.

      INTEGER ::                                 LAMBDA_ITER

      LOGICAL, DIMENSION(:), ALLOCATABLE ::      NOFILE

      logical ::                                 each_ali = .false.

      LOGICAL ::                                 VEL_FIELD_FROM_FILE

      REAL*8, DIMENSION(:), ALLOCATABLE ::       HEIGHT

      integer ::                                 NDPMIN ! THE DEPTH POINT OF TEMPERATURE MINIMUM

      LOGICAL ::                                 OLDSTART

      INTEGER ::                                 NLINE

      real*8, allocatable ::                     ABXYZ(:), ABXYZn(:, :)

      integer ::                                 negintl

!     nlte marking variables

      logical, dimension(:), allocatable ::      nlte_lev, nlte_ele

!      integer ::                                 N_nlte, lastind_nlte, natom_nlte

      integer, dimension(:), allocatable ::      idx_orig

      real*8,                allocatable ::      ABXYZ_nlte(:), ABXYZn_nlte(:, :)

      logical ::                                 LBKG           ! KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTION (ODF)
      integer ::                                 XLBKG1, XLBKG2 ! WAVELENTH RANGE FOR THE ODF

!     FIOSS VARIABLES

      INTEGER ::                                 NUMTRA

      INTEGER ::                                 NTC

      INTEGER, ALLOCATABLE, DIMENSION(:) ::      NTI

      INTEGER, ALLOCATABLE, DIMENSION(:) ::      ELID, CHARGE, SWL, SWU

      REAL*8,  ALLOCATABLE, DIMENSION(:) ::      WAVTRA, AUL

      REAL*8,  ALLOCATABLE, DIMENSION(:, :) ::   NTPL, NTPU, NTDL, NTDU

!     ------------------------------------------------------------------------------------
!     variables and parameters for the odf_table module

      integer :: ntemp, nrhox

      real*8, allocatable, dimension(:) :: tabt, tabp

      integer, parameter :: nsizebig = 90
      integer, parameter :: isubbin =  10

      integer(kind = 2), allocatable, dimension(:, :, :, :) :: iodfrecvb
      real*8,            allocatable, dimension(:)          :: freqgrid

      end module
