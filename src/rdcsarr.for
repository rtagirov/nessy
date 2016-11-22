      module MOD_RDCSARR

      contains

      SUBROUTINE RDCSARR(LEVEL,J,N,WAVARR,SIGARR,NDIM,NFDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)
CMH   CALLED BY DATOM
CMH   READS WAVENUMBER AND CROSS SECTIONS FOR EACH EXPLICIT LEVEL INTO AN ARRAY
CMH   NDIM : NUMBER OF LEVELS
CMH   NFDIM: NUMBER OF FREQUENCIES
CMH   WAVARR: WAVENUMBERS FOR EACH LEVEL
CMH   SIGARR: CROSS SECTIONS FOR EACH LEVEL
      real*8,intent(inout)  :: WAVARR,SIGARR
      integer,intent(in) :: J,N,NDIM,NFDIM
      CHARACTER*10,intent(in) :: LEVEL(NDIM)

      INTEGER IOSTATUS
      CHARACTER*10,FILENAME
      real*8 WLTH, SIGMA
      k=1
      SIGARR(J,:)=0.
      WAVARR(J,:)=0.
      FILENAME = LEVEL(J)
      open(unit=1, file=FILENAME,STATUS='OLD',IOSTAT=IOSTATUS,err=888,
     * action='read')
      do  while (IOSTATUS .eq. 0)
        read(unit=1, fmt=*, IOSTAT=IOSTATUS), WLTH, SIGMA
        if (wlth .eq. 0.0) then
            print *,'RDCSARR: PROBLEM CROSS SECTION INPUT!'
            PAUSE
          endif
          wavarr(J,K)=1.d8/wlth
          sigarr(J,K)=sigma
          k=k+1
      end do
888   continue
      if (k .eq. 1) then
        print *,'RDCSARR: PROBLEM with ',FILENAME
        pause
      endif
      close (unit=1)
c     endif
    1 continue
      return
      END subroutine
      end module
