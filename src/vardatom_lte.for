      module vardatom_lte

      integer                                    :: natom_lte, lis_num

      integer,      dimension(:),    allocatable :: eleatnum_lte, eleisnum_lte

      integer,      dimension(:),    allocatable :: lis_anum, lis_cnum, lis_lnum

      integer,      dimension(:, :), allocatable :: lis_weight

      real*8,       dimension(:, :), allocatable :: lis_levien

      character*10, dimension(:),    allocatable :: lis_name

      real*8,       dimension(:),    allocatable :: abxyz_lte

      real*8,       dimension(:, :), allocatable :: abxyzn_lte

      end module
