      module odf_table

      contains

      subroutine read_odf_table()

!     reads the ODF table as given by DFSYNTHE calculations
!     the ODF is given for a fixed value of turbulent velocity

      use common_block
!      use phys

      implicit none

      integer :: j, inu, istep, it, ip

      real*8, allocatable, dimension(:) :: temp, pres

!      real*8 :: rho ! density for a given pair of pressure and temperature

      open(unit = 1409, file = 'odf.table.grid')

      read(1409, *) numt
      read(1409, *) nump

      allocate(temp(numt), pres(nump))
      allocate(tabt(numt), tabp(nump))

      read (1409, *) (temp(j), j = 1, numt)
      read (1409, *) (pres(j), j = 1, nump)

      close(unit = 1409)

      tabt = dlog10(temp)
      tabp = dlog10(pres)

      allocate(odf(nsubbins, nbins, nump, numt))
      allocate(wvlgrid(nbins + 1))

      open(unit = 1408, file = 'odf.table')

      do inu = 1, nbins

         read(1408, *) wvlgrid(inu), wvlgrid(inu + 1)

            do ip = 1, nump

               do it = 1, numt

!                  rho = pres(ip) * apm / boltz / temp(it)
!                  rho = pres(ip) * apm / (boltz * temp(it) + apm * 1d+10 / 2.0d0)

!                  print*, 'lalala', apm

!                  stop

                  read(1408, *) (odf(istep, inu, ip, it), istep = 1, nsubbins)

!                  odf(1 : nsubbins, inu, ip, it) =  1d+3 * dlog10(rho * 1.0d+1**(odf(1 : nsubbins, inu, ip, it) / 1d+3))

               enddo

            enddo

      enddo

      close(unit = 1408)

      deallocate(temp)
      deallocate(pres)

      return
 
      end subroutine


      subroutine odf_interpolation_coef(n, ne, T, idxt, idxp, c1, c2, c3, c4)

      use phys
      use common_block

      implicit none

!--------------------------------- CONSTANTS ---------------------------

!     tenlog is the conversion from log10 to ln
      real*8, parameter :: tenlog = 2.30258509299405d0
      real*8, parameter :: ttenlg = 0.001d0 * tenlog

!------------------------------- IN-OUT -------------------------------

      real*8, intent(in), dimension(dpn) ::   T, n, ne

      integer, intent(out), dimension(dpn) :: idxt, idxp

      real*8, intent(out), dimension(dpn) ::  c1, c2, c3, c4

!------------------------------- LOCAL VARIABLES -----------------------

      real*8  :: p, plog, tlog, x, y

      integer :: j, ip, it

!---------------------------------- EXECUTION --------------------------

!      if (.not. allocated(idxt)) allocate(idxt(dpn))
!      if (.not. allocated(idxp)) allocate(idxp(dpn))

!      if (.not. allocated(c1)) allocate(c1(dpn))
!      if (.not. allocated(c2)) allocate(c2(dpn))
!      if (.not. allocated(c3)) allocate(c3(dpn))
!      if (.not. allocated(c4)) allocate(c4(dpn))

      do j = 1, dpn

         tlog = min(max(log10(T(j)), tabt(1)), tabt(numt))

         it = 2

         do while (tabt(it) .le. tlog .and. it .lt. numt)

            it = it + 1

         enddo

         p = (n(j) + ne(j)) * boltz * T(j)

         plog = min(max(log10(p), tabp(1)), tabp(nump))

         ip = 2

         do while (tabp(ip) .le. plog .and. ip .lt. nump)

            ip = ip + 1

         enddo

         idxp(j) = ip
         idxt(j) = it

         x = (tlog - tabt(it - 1)) / (tabt(it) - tabt(it - 1))
         y = (plog - tabp(ip - 1)) / (tabp(ip) - tabp(ip - 1))

!        the coefficients are scaled back by 0.001 by ttenlg because opacities 
!        are read from odf.table as 1000 * log10(opacity)
         c1(j) = (1.0d0 - x) * (1.0d0 - y) * ttenlg
         c2(j) = (1.0d0 - x) * y * ttenlg
         c3(j) = x * (1.0d0 - y) * ttenlg
         c4(j) = x * y * ttenlg

      enddo

      return

      end subroutine


      subroutine odf_interpolation(xlam, n, idxt, idxp, c1, c2, c3, c4, linop)

      use file_operations
      use common_block
      use utils
      use phys

      implicit none

!     ODF interpolation for all depth points and single wavelength value (xlam)
!     using the precalculated interpolation coefficients (see above, odf_interpolation_coef)
!     ODF has been read into the odf array from the file odf.table (see above subroutine read_odf_table)
!     The odf array is a global variable (see comblock.for)

!------------------------------- IN-OUT -------------------------------

      real*8, intent(in)                  :: xlam

      real*8, intent(in), dimension(dpn)  :: n

      integer, intent(in), dimension(dpn) :: idxt, idxp

      real*8, intent(in), dimension(dpn)  :: c1, c2, c3, c4

      real*8, intent(out), dimension(dpn) :: linop

!------------------------------- LOCAL VARIABLES -----------------------

      real*8, dimension(nsubbins + 1) :: subgrid

      real*8, dimension(dpn) :: rho

      real*8 :: delta

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

         it = idxt(j)
         ip = idxp(j)

         linop(j) = dexp(c1(j) * dble(odf(sbn, bn, ip - 1, it - 1)) +
     &                   c2(j) * dble(odf(sbn, bn, ip  ,   it - 1)) +
     &                   c3(j) * dble(odf(sbn, bn, ip - 1, it)) +
     &                   c4(j) * dble(odf(sbn, bn, ip  ,   it)))

      enddo

!     n is the total particle concentration (everything except electrons)
!     apm is the average particle mass
!     rho is density
      rho = n * apm

      linop = rho * linop

!      call open_to_append(1435, 'linop.out')

!      do j = 1, dpn

!         write(1435, '(e15.7,3(1x,i4),1x,e15.7)') xlam, bn, sbn, j, linop(j)

!      enddo

!      close(1435)

      linop(1 : ndpmin) = linop(ndpmin)
!      linop(1 : ndpmin) = 0.0d0
!      linop = linop * 1.0d-10

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

      if (n == nbins + 1)    cond = grid(idx + 1) .le. w
      if (n == nsubbins + 1) cond = grid(idx + 1) .lt. w

      do while (cond)

         idx = idx + 1

         if (n == nbins + 1)    cond = grid(idx + 1) .le. w
         if (n == nsubbins + 1) cond = grid(idx + 1) .lt. w

      enddo

      return

      end function

      end module
