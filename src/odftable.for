      module odf_table

      contains

      subroutine read_odf_table()

      use common_block

      implicit none

      integer :: j, inu, istep, it, ip

      open(unit = 1409, file = 'odf.table.grid')

      read(1409, *) numt
      read(1409, *) nump

      allocate(tabt(numt), tabp(nump))

      read (1409, *) (tabt(j), j = 1, numt)
      read (1409, *) (tabp(j), j = 1, nump)

      close(unit = 1409)

      tabt = log10(tabt)
      tabp = log10(tabp)

      allocate(odf(nsubbins, nbins, nump, numt))
      allocate(wvlgrid(nbins + 1))

      open(unit = 1408, file = 'odf.table')

      do inu = 1, nbins

         read(1408, *) wvlgrid(inu), wvlgrid(inu + 1)

            do ip = 1, nump

               do it = 1, numt

                  read(1408, *) (odf(istep, inu, ip, it), istep = 1, nsubbins)

               enddo

            enddo

      enddo

      close(unit = 1408)

      return
 
      end subroutine


      subroutine odf_interpolation_coef(entot, T)

      use phys
      use common_block

      implicit none

!--------------------------------- CONSTANTS ---------------------------

!     tenlog is the conversion from log10 to ln
      real*8, parameter :: tenlog = 2.30258509299405d0
      real*8, parameter :: ttenlg = 0.001d0 * tenlog

!------------------------------- IN-OUT -------------------------------

      real*8, intent(in), dimension(dpn) :: T, entot

!------------------------------- LOCAL VARIABLES -----------------------

      real*8  :: p, plog, tlog, x, y

      integer :: j, ip, it

!---------------------------------- EXECUTION --------------------------

      allocate(idx_temp(dpn), idx_pres(dpn))
      allocate(co1(dpn), co2(dpn), co3(dpn), co4(dpn))

      do j = 1, dpn

         tlog = min(max(log10(T(j)), tabt(1)), tabt(numt))

         it = 2

         do while (tabt(it) .le. tlog .and. it .lt. numt)

            it = it + 1

         enddo

         p = entot(j) * boltz * T(j)

         plog = min(max(log10(p), tabp(1)), tabp(nump))

         ip = 2

         do while (tabp(ip) .le. plog .and. ip .lt. nump)

            ip = ip + 1

         enddo

         idx_pres(j) = ip
         idx_temp(j) = it

         x      = (tlog - tabt(it - 1)) / (tabt(it) - tabt(it - 1))
         y      = (plog - tabp(ip - 1)) / (tabp(ip) - tabp(ip - 1))

!        the coefficients are scaled back by 0.001 by ttenlg because opacities 
!        are read from odf.table as 1000 * log10(opacity)
         co1(j) = (1.0d0 - x) * (1.0d0 - y) * ttenlg
         co2(j) = (1.0d0 - x) * y * ttenlg
         co3(j) = x * (1.0d0 - y) * ttenlg
         co4(j) = x * y * ttenlg

      enddo

      return

      end subroutine


      subroutine odf_interpolation(xlam, linop)

      use common_block
      use utils

      implicit none

!     Assumes that vturb is constant and that the opacity file is given only for that vturb
!     ODF have been already read into odf array (which is a global variable, see comblock.for) 
!     from the odf.table file

!------------------------------- IN-OUT -------------------------------

      real*8, intent(in)                   :: xlam

      real*8, intent(out), dimension(dpn)  :: linop

!------------------------------- LOCAL VARIABLES -----------------------

      real*8, dimension(nsubbins + 1) :: subgrid

      real*8  :: delta

      integer :: k, j, ip, it, bn, sbn

!---------------------------------- EXECUTION --------------------------

      linop = 0.0d0

      if (xlam < xlbkg1) return
      if (xlam > xlbkg2) return

!     bn is the bin number
      bn = bin_index(nbins + 1, wvlgrid, xlam / 10.0)

      delta = (wvlgrid(bn + 1) - wvlgrid(bn)) / nsubbins

      do k = 1, nsubbins + 1

         subgrid(k) = wvlgrid(bn) + (k - 1) * delta

      enddo

      sbn = bin_index(nsubbins + 1, subgrid, xlam / 10.0)

      do j = 1, dpn

         it = idx_temp(j)
         ip = idx_pres(j)

         linop(j) = dexp(co1(j) * dble(odf(sbn, bn, ip - 1, it - 1)) +
     &                   co2(j) * dble(odf(sbn, bn, ip  ,   it - 1)) +
     &                   co3(j) * dble(odf(sbn, bn, ip - 1, it)) +
     &                   co4(j) * dble(odf(sbn, bn, ip  ,   it)))

      enddo

      linop(1 : ndpmin) = linop(ndpmin)

      if (any(linop < 0.0d0)) then

          print*, xlam, 'linop = ', linop

          call error('odf_interpolation: negative opacity')

      endif

      return

      end subroutine


      function bin_index(n, grid, w) result(idx)

      use common_block

      implicit none

      integer              :: n

      real*8, dimension(n) :: grid

      real*8               :: w

      integer              :: idx

      logical              :: cond

      idx = 1

!     grid contains borders of the bins, such that grid(idx + 1) corresponds to the bin's upper wavelength
!     for wvlgrid this is read from the odf.table, which gives the beginning and the end wavelength of each bin
!     for subgrid this is calculated in the odf_interpolation subroutine

      if (n == nbins)    cond = grid(idx + 1) .le. w
      if (n == nsubbins) cond = grid(idx + 1) .lt. w

      do while (cond)

         idx = idx + 1

         if (n == nbins)    cond = grid(idx + 1) .le. w
         if (n == nsubbins) cond = grid(idx + 1) .lt. w

      enddo

      return

      end function

      end module
