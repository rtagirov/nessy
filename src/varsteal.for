      module varsteal

!     variables shared between steal, como and etl

      real*8,  allocatable, dimension(:)    ::  phi, pweight

      integer, allocatable, dimension(:)    ::  levelpl, itne, iwarn

      real*8,  allocatable, dimension(:)    ::  opa, eta

      real*8,  allocatable, dimension(:)    ::  thomson

      real*8,  allocatable, dimension(:, :) ::  depart

      end module
