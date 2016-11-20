      function nan()
      real*8 nan
      real*8,parameter :: inf=1/0d0
      real*8,parameter :: nand=inf/inf
      nan=nand
      end
