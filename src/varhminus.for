      module varhminus

      use params_array

!     variables shared between wrstart, steal, wrcont, como and etl

      real*8, dimension(NFDIM)                :: XLAMBDA, AKEY, FWEIGHT
      
      real*8, allocatable, dimension(:)       :: radius
      real*8, allocatable, dimension(:)       :: entot
      real*8, allocatable, dimension(:)       :: T

      real*8, allocatable, dimension(:)       :: VELO, GRADI

      real*8, allocatable, dimension(:)       :: TAUROSS

      real*8, allocatable, dimension(:)       :: EMFLUX
      real*8, allocatable, dimension(:)       :: HTOT, GTOT, XTOT, ETOT

      real*8, allocatable, dimension(:)       :: P, Z

      real*8, allocatable, dimension(:)       :: XJC
      real*8, allocatable, dimension(:, :)    :: XJCARR, XJL
      real*8, allocatable, dimension(:, :)    :: EDDI
      real*8, allocatable, dimension(:, :, :) :: EDDARR

      real*8, allocatable, dimension(:)       :: ENLTE, RNE
      real*8, allocatable, dimension(:, :)    :: POPNUM, POP1, POP2, POP3
      real*8, allocatable, dimension(:, :)    :: WCHARM

      end module
