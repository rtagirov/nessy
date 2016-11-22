      module MOD_RDCSARR

      contains

      SUBROUTINE RDCSARR(LEVEL, J, WAVARR, SIGARR, N, NFDIM)

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION WAVARR(N, NFDIM), SIGARR(N, NFDIM)
CMH   CALLED BY DATOM
CMH   READS WAVENUMBER AND CROSS SECTIONS FOR EACH EXPLICIT LEVEL INTO AN ARRAY
CMH   N :     NUMBER OF LEVELS
CMH   NFDIM:  MAX NUMBER OF FREQUENCIES
CMH   WAVARR: WAVENUMBERS FOR EACH LEVEL
CMH   SIGARR: CROSS SECTIONS FOR EACH LEVEL
      integer, intent(in) :: J, N, NFDIM
      CHARACTER*10, intent(in) :: LEVEL(N)

      real*8, intent(out) :: WAVARR, SIGARR

      INTEGER IOSTATUS
      CHARACTER*10, FILENAME
      real*8 WLTH, SIGMA

      k = 1

      SIGARR(J, :) = 0.
      WAVARR(J, :) = 0.

      FILENAME = LEVEL(J)

      open(unit = 1, file=FILENAME, STATUS='OLD', IOSTAT=IOSTATUS, err=888, action='read')

      do while (IOSTATUS .eq. 0)

         read(unit = 1, fmt=*, IOSTAT=IOSTATUS), WLTH, SIGMA

         if (wlth .eq. 0.0) then

             print *,'RDCSARR: PROBLEM CROSS SECTION INPUT!'

             stop

         endif

!         print*, 'rdcsarr: ', level(j), wlth, sigma

         wavarr(J, K) = 1.0D8 / wlth
         sigarr(J, K) = sigma

         k = k + 1

      enddo

888   continue

      if (k .eq. 1) then

        print *,'RDCSARR: PROBLEM with ',FILENAME

        stop

      endif

      close(unit = 1)

    1 continue

      return

      END subroutine

      end module
