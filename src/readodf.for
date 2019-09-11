module readodf

contains

subroutine read_odf(nb, nw, np, nT, odf)

use common_block

implicit none

integer, intent(in):: nb, nw, np, nT

real*8, intent(out), allocatable, dimension(:, :, :) :: odf

real*8 :: x

integer :: i, j, k

print*, 'check start'

allocate(odf(nw, np, nT))

open(unit = 1408, file = 'atlas9.odf')

do k = 1, nb

    read(nl, *) x, x

    do i = 1, np

       do j = 1, nT

            read(nl, *) odf(10 * (k - 1) + 1 : 10 * k, i, j)

       enddo 

    enddo

enddo

close(unit = 1408)

odf = 10.0**(odf / 1000.0)

print*, 'check end'

return

end subroutine

end module
