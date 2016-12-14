      module varsteal

      use params_array

!     variables shared between steal, como and etl

      real*8,  dimension(NFLDIM)            ::  phi, pweight

      integer, allocatable, dimension(:)    ::  levelpl, nfedge, itne, iwarn

      real*8,  allocatable, dimension(:)    ::  en

      real*8,  allocatable, dimension(:)    ::  opa, eta

      real*8,  allocatable, dimension(:)    ::  thomson, tauthom

      real*8,  allocatable, dimension(:)    ::  opac, etac

      real*8,  allocatable, dimension(:)    ::  dopa, deta

      real*8,  allocatable, dimension(:)    ::  expfac

      real*8,  allocatable, dimension(:, :) ::  sigmaki, depart

      end module
