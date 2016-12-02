      module varhminus

      use params_array

      real*8, dimension(NFDIM) :: XLAMBDA, AKEY, FWEIGHT
      
      real*8, allocatable, dimension(:) :: radius
      real*8, allocatable, dimension(:) :: entot
      real*8, allocatable, dimension(:) :: T

      real*8, allocatable, dimension(:) :: XJC
      real*8, allocatable, dimension(:) :: TAUROSS
      real*8, allocatable, dimension(:) :: RNE
      real*8, allocatable, dimension(:) :: VELO
      real*8, allocatable, dimension(:) :: GRADI
      real*8, allocatable, dimension(:) :: EMFLUX
      real*8, allocatable, dimension(:) :: ENLTE
      real*8, allocatable, dimension(:) :: P
      real*8, allocatable, dimension(:) :: Z
      real*8, allocatable, dimension(:) :: HTOT, GTOT, XTOT, ETOT

      real*8, allocatable, dimension(:, :) :: XJCARR
      real*8, allocatable, dimension(:, :) :: XJL
      real*8, allocatable, dimension(:, :) :: EDDI
      real*8, allocatable, dimension(:, :) :: POPNUM, POP1, POP2, POP3
      real*8, allocatable, dimension(:, :) :: WCHARM

      real*8, allocatable, dimension(:, :, :) :: EDDARR

      end module
