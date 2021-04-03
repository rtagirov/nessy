      module mod_synopa

      contains

      subroutine synopa(WAVARR, SIGARR, N, NFDIM)

      use MOD_SYNSUBM
      use SYNTHP_CONT
      use OPINT, only: OPATOT, ETATOT
      
      use auxfioss
      use mod_opac

!     interface opacity routine

      implicit none

      integer, intent(in) :: N, NFDIM
      real*8,  intent(in) :: WAVARR, SIGARR
      real*8 ::              ABSO, EMIS
    
      integer :: ID, IJ, i

      real*8 :: opac_start, opac_finish

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/MODELP.FOR'

      DIMENSION WAVARR(N, NFDIM), SIGARR(N, NFDIM)

CMH   MARGIT HABERREITER
CMH   IJ:    INDEX OF FREQUENCT
CMH   NFREQ: NUMBER OF FREQUENCY POINTS
CMH   INPUT PARAMETERS:

CMH   OUTPUT PARAMETERS:
CMH   OPATOT: NEW ARRAY FOR OPACITY for each frequency nfreq and depth point nd
CMH   ETATOT: NEW ARRAY FOR EMISSIVITY for each frequency nfreq and depth point nd
      dimension abso(mfreq), emis(mfreq)

!      call system('rm -vf linop_loops.time')

!      open(unit = 18765, file = 'linop_loops.time')

!     see linop.for; look for 18765
!      write(18765, '(7(3x,A),3(6x,A))'), 'dp', 'total', 'ommited', 'diff',
!     $                                   'voigt', 'svc', 'dvc', 
!     $                                   'cycle1', 'cycle2', 'cycle3'

      call cpu_time(opac_start)

      do 20 id = 1, nd

         call opac(id, 1, abso, emis, WAVARR, SIGARR, N, NFDIM)

!        loop over frequency points
         do ij = 1, nfreq

cws May-24-96: Ivan's routine calculated low to high freqency
c         wlam(nvopa-kopa+1) = xlam * (1. - vel*clkm)
c         freq(nvopa-kopa+1) = cl8/wlam(nvopa-kopa+1)
            opatot(nfreq-ij+1,id)=abso(ij)
            etatot(nfreq-ij+1,id)=emis(ij)

!            print*, 'synopa test', abso(ij)

            if (id.le.1) etatot(nfreq-ij+1,id)=emisc(getContIdx(ij))

         enddo

   20 continue

      call cpu_time(opac_finish)

!      close(18765)

      print*, 'synopa: opac execution time = ', (opac_finish - opac_start)

  310 format (1pe12.5, 1pe12.5)

      return

      end subroutine

      end module
