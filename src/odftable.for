      module odf_table

      contains

!      subroutine read_TP_grid(ntemp, nrhox, tabt, tabp)
      subroutine read_tp_grid()

      use common_block

      implicit none

!      integer, intent(out) :: ntemp, nrhox

!      real*8, allocatable, dimension(:) :: tabt, tabp

      integer :: j

      open(unit = 1409, file = 'tp.grid')

      read(1409, *) ntemp
      read(1409, *) nrhox

      allocate(tabt(ntemp), tabp(nrhox))

      read (1409, *) (tabt(j), j = 1, ntemp)
      read (1409, *) (tabp(j), j = 1, nrhox)

      close(unit = 1409)

      return

      end subroutine


!      subroutine read_odf_table(nsizebig, isubbin, ntemp, nrhox, freqgrid, iodfrecvb)
      subroutine read_odf_table()

      use common_block

      implicit none

!      integer, intent(in) :: nsizebig, isubbin, nrhox, ntemp

!      integer(kind = 2), intent(out), allocatable, dimension(:, :, :) :: iodfrecvb
!      real*8,            intent(out), allocatable, dimension(:)       :: freqgrid

      integer :: inu, istep, it, ip

!      print*, 'check start'

      allocate(iodfrecvb(isubbin, nsizebig, nrhox, ntemp))
      allocate(freqgrid(nsizebig + 1))

      open(unit = 1408, file = 'odf.table')

      do inu = 1, nsizebig

         read(1408, *) freqgrid(inu), freqgrid(inu + 1)

            do ip = 1, nrhox

               do it = 1, ntemp

                  read(1408, *) (iodfrecvb(istep, inu, ip, it), istep = 1, isubbin)

               enddo

            enddo

      enddo

      close(unit = 1408)

!      print*, 'check end'

!      stop

      return
 
      end subroutine


      subroutine interPT(n, stepwt)

      implicit none

!     Assumes that vturb is constant and that the opacity file is given only for that vturb
!     ODF have been already read into kapw which is a global variable from the odf.table file

!---------------------------- DUMMY VARIABLES --------------------------

      integer, intent(in) :: n
      integer             :: n1

!--------------------------------- CONSTANTS ---------------------------



      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!  tenlog is conversion from log10 to ln ---   parameter ( tenlog = 2.30258509299405d0)
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision, intent(out) :: stepwt
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd),
     &                  plog, tl, wave, x, y
      integer    j, ip, ipj(maxd), it, itemp1, itj(maxd), iwl
      save       co1, co2, co3, co4,  ipj, itj, iwl
!
!----------------------------- INITIALIZATION --------------------------
!
      data itemp1 / 0 /
!
!---------------------------------- EXECUTION --------------------------
!
! - iwl is the loop index over all ibin s
!
! isubbin - nsteps
! wlend it the freq where the bin finishes, so inifreset(2:ibin+1)
! global binwith variable, sbwt,  and globale boarders of subbins sbvalues
     if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
! tlog(j) -- log(T(j)) where j- is the depth point of your atmosphere !
!
         do j = 1, nrhox
            tl = min(max(tlog(j) / tenlog, tabt(1)), tabt(numt) )
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. numt)
               it = it + 1
            end do
!
            plog = min(max(log10(p(j)), tabp(1)), tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. numpres)
               ip = ip + 1
            end do
!
            ipj(j) = ip
            itj(j) = it
            x      = (tl   - tabt(it-1)) / (tabt(it) - tabt(it-1))
            y      = (plog - tabp(ip-1)) / (tabp(ip) - tabp(ip-1))
            co1(j) = (1.0d0 - x) * (1.0d0 - y) * ttenlg
            co2(j) = (1.0d0 - x) * y * ttenlg
            co3(j) = x * (1.0d0 - y) * ttenlg
            co4(j) = x * y * ttenlg
!
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
! inifreset contains bin boarders, such that iwl+1 corresponds to upper bin wavelength 
! thes can be read from the ODF file, as for each bin it gives the beginning and end wavelengt
!
      do while (inifreset(iwl+1) .le. wave)
         iwl = iwl + 1
      end do
!    we start with the smallest sub-bin
      n1= isubbin+1 -n
      stepwt = sbwt(n1,iwl)
! but now there are globally known!  (sbwt(n, iwl) = wt(n), nsteps = isubbin)
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(kapw(n1, iwl, ip-1, it-1)) +
     &             co2(j) * dble(kapw(n1, iwl, ip  , it-1)) +
     &             co3(j) * dble(kapw(n1, iwl, ip-1, it)) +
     &             co4(j) * dble(kapw(n1, iwl, ip  , it)) )
         alines(j) = a
      end do

      return

      end subroutine

      end module
