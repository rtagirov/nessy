      module vardatom_full

!     datom_full file variables (both for fioss and hminus)

      integer                                   :: N, NATOM, LASTIND

      integer, allocatable, dimension(:)        :: KODAT, NOM, NFIRST, NLAST, INDLOW, INDNUP, NCHARG, MAINQN

      integer, allocatable, dimension(:)        :: eleatnum, levatnum

      real*8,  allocatable, dimension(:)        :: ATMASS, ALPHA, EION, ELEVEL, SEXPO, WEIGHT, STAGE

      real*8,  allocatable, dimension(:, :)     :: ALTESUM, EINST, SIGARR, WAVARR

      real*8,  allocatable, dimension(:, :, :)  :: COCO

      character*2, allocatable, dimension(:)    :: symbol

      character*4, allocatable, dimension(:, :) :: keycol

      character*8, allocatable, dimension(:)    :: agaunt

      character*10, allocatable, dimension(:)   :: level, element

      end module
