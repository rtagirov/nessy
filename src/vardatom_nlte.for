      module vardatom_nlte

!     variables for subroutine datom (both for fioss and hminus); mode = 'nlte'

      integer                                   :: N_nlte, natom_nlte, lastind_nlte

      integer, allocatable, dimension(:)        :: kodat_nlte, nom_nlte, nfirst_nlte, nlast_nlte,
     $                                             indlow_nlte, indnup_nlte, ncharg_nlte, mainqn_nlte

      integer, allocatable, dimension(:)        :: eleatnum_nlte, levatnum_nlte

      real*8,  allocatable, dimension(:)        :: atmass_nlte, alpha_nlte, eion_nlte, elevel_nlte,
     $                                             sexpo_nlte, weight_nlte, stage_nlte

      real*8,  allocatable, dimension(:, :)     :: altesum_nlte, einst_nlte, sigarr_nlte, wavarr_nlte

      real*8,  allocatable, dimension(:, :, :)  :: coco_nlte

      character*2, allocatable, dimension(:)    :: symbol_nlte

      character*4, allocatable, dimension(:, :) :: keycol_nlte

      character*8, allocatable, dimension(:)    :: agaunt_nlte

      character*10, allocatable, dimension(:)   :: level_nlte, element_nlte

      end module
