      module varhminus

!     variables shared between wrstart, steal, wrcont, como and etl

      integer, parameter                      :: NFDIM =   2000

      real*8, dimension(NFDIM)                :: XLAMBDA, AKEY, FWEIGHT
      
      real*8, allocatable, dimension(:)       :: radius
      real*8, allocatable, dimension(:)       :: entot
      real*8, allocatable, dimension(:)       :: T

      real*8, allocatable, dimension(:)       :: VELO, GRADI

      real*8, allocatable, dimension(:)       :: TAUROSS

      real*8, allocatable, dimension(:)       :: EMFLUX
      real*8, allocatable, dimension(:)       :: HTOT, GTOT, XTOT, ETOT

      real*8, allocatable, dimension(:)       :: P, Z

      real*8, allocatable, dimension(:)       :: xjc
      real*8, allocatable, dimension(:, :)    :: xjc2, XJL
      real*8, allocatable, dimension(:, :)    :: EDDI
      real*8, allocatable, dimension(:, :, :) :: EDDARR

      real*8, allocatable, dimension(:)       :: RNE
      real*8, allocatable, dimension(:, :)    :: POPNUM, POP1, POP2, POP3
      real*8, allocatable, dimension(:, :)    :: WCHARM

      real*8, allocatable, dimension(:, :)    :: scafac, absfac

      real*8, allocatable, dimension(:, :)    :: U, VJL

      real*8, allocatable, dimension(:, :)    :: sigmaki

      real*8, allocatable, dimension(:)       :: VL, HNU

      character*10, allocatable, dimension(:) :: mainpro, mainlev

      end module
