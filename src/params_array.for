      module PARAMS_ARRAY

      integer, parameter :: MAXATOM = 90  ! DEFINE ARRAY DIMENSIONS
      integer, parameter :: NDIM =    301 ! NUMBER OF LEVELS
      integer, parameter :: NDIMP2 =  NDIM + 2
      integer, parameter :: NFDIM =   2000
      integer, parameter :: NFLDIM =  79

      integer, parameter :: NDDIM =  200  ! DEPTH POINTS
      integer, parameter :: NPDIM =  208  ! DEPTH POINTS + 8 ADDITIONAL POINTS
      integer, parameter :: MAXIND = 1000 

      end module

! Achtung: fioss.for definiert seine eigenen Parameter
