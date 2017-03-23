      module vardatom

!     DATOM file variables (both for fioss and hminus)

      integer, allocatable, dimension(:)        :: KODAT, NOM, NFIRST, NLAST, INDLOW, INDNUP, NCHARG, MAINQN

      real*8,  allocatable, dimension(:)        :: ATMASS, ALPHA, EION, ELEVEL, SEXPO, WEIGHT, STAGE

      real*8,  allocatable, dimension(:, :)     :: ALTESUM, EINST, SIGARR, WAVARR

      real*8,  allocatable, dimension(:, :, :)  :: COCO

      character*2, allocatable, dimension(:)    :: symbol

      character*4, allocatable, dimension(:, :) :: keycol

      character*8, allocatable, dimension(:)    :: agaunt

      character*10, allocatable, dimension(:)   :: level, element

      integer, allocatable, dimension(:)        :: kodat_nlte, nom_nlte, nfirst_nlte, nlast_nlte,
     $                                             indlow_nlte, indnup_nlte, ncharg_nlte, mainqn_nlte

      real*8,  allocatable, dimension(:)        :: atmass_nlte, alpha_nlte, eion_nlte, elevel_nlte,
     $                                             sexpo_nlte, weight_nlte, stage_nlte

      real*8,  allocatable, dimension(:, :)     :: altesum_nlte, einst_nlte, sigarr_nlte, wavarr_nlte

      real*8,  allocatable, dimension(:, :, :)  :: coco_nlte

      character*2, allocatable, dimension(:)    :: symbol_nlte

      character*4, allocatable, dimension(:, :) :: keycol_nlte

      character*8, allocatable, dimension(:)    :: agaunt_nlte

      character*10, allocatable, dimension(:)   :: level_nlte, element_nlte

      end module vardatom
