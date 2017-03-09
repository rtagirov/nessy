      module mod_synopa

      contains

      subroutine synopa(WAVARR, SIGARR, N, NF)

      use MOD_SYNSUBM
      use SYNTHP_CONT
      use OPINT, only: OPATOT, ETATOT
      
      use auxfioss
      use mod_opac

!     interface opacity routine

      implicit none

      integer, intent(in) :: N, NF
      real*8,  intent(in) :: WAVARR, SIGARR
      real*8 ::              ABSO, EMIS
    
      integer :: ID, IJ, i

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/MODELP.FOR'

      DIMENSION WAVARR(N, NF), SIGARR(N, NF)

CMH   MARGIT HABERREITER
CMH   IJ:    INDEX OF FREQUENCT
CMH   NFREQ: NUMBER OF FREQUENCY POINTS
CMH   INPUT PARAMETERS:

CMH   OUTPUT PARAMETERS:
CMH   OPATOT: NEW ARRAY FOR OPACITY for each frequency nfreq and depth point nd
CMH   ETATOT: NEW ARRAY FOR EMISSIVITY for each frequency nfreq and depth point nd
      dimension abso(mfreq), emis(mfreq)

      do 20 id = 1, nd

         call opac(id, 1, abso, emis, WAVARR, SIGARR, N, NF)

!        loop over frequency points
         do ij = 1, nfreq

            print*, 'freq = ', freq(ij)

cws May-24-96: Ivan's routine calculated low to high freqency
c         it is (in QUANT): 
c         wlam(nvopa-kopa+1) = xlam * (1. - vel*clkm)
c         freq(nvopa-kopa+1) = cl8/wlam(nvopa-kopa+1)
            opatot(nfreq-ij+1,id)=abso(ij)
            etatot(nfreq-ij+1,id)=emis(ij)

            if (id.le.1) etatot(nfreq-ij+1,id)=emisc(getContIdx(ij))

         enddo

   20 continue

  310 format (1pe12.5, 1pe12.5)

      return

      end subroutine

      end module
